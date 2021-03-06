package main

import (
	"flag"
	"fmt"
	"github.com/kortschak/BioGo/align/pals"
	"github.com/kortschak/BioGo/align/pals/filter"
	"github.com/kortschak/BioGo/morass"
	"github.com/kortschak/BioGo/seq"
	"io"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"runtime/pprof"
	"strconv"
	"strings"
	"time"
)

const (
	timeFormat = "20060102150405-Mon"
)

var (
	pid           = os.Getpid()
	mem           *uintptr
	profile       *os.File
	queryName     string
	targetName    string
	selfCompare   bool
	sameStrand    bool
	outFile       string
	maxK          int
	minHitLen     int
	minId         float64
	dpMinHitLen   int
	dpMinId       float64
	tubeOffset    int
	tmpDir        string
	tmpChunk      int
	tmpConcurrent bool
	threads       int
	maxMem        uint64
	logToFile     bool
	debug         bool
	verbose       bool
	cpuprofile    string
	logger        *log.Logger
)

func init() {
	flag.StringVar(&queryName, "query", "", "Filename for query sequence.")
	flag.StringVar(&targetName, "target", "", "Filename for target sequence.")
	flag.BoolVar(&selfCompare, "self", false, "Is this a self comparison?")
	flag.BoolVar(&sameStrand, "same", false, "Only compare same strand")

	flag.StringVar(&outFile, "out", "", "File to send output to.")

	flag.IntVar(&maxK, "k", -1, "Maximum kmer length (negative indicates automatic detection based on architecture).")
	flag.IntVar(&minHitLen, "filtlen", 400, "Minimum hit length for filter.")
	flag.Float64Var(&minId, "filtid", 0.94, "Minimum hit identity for filter.")
	flag.IntVar(&dpMinHitLen, "dplen", 0, "Minimum hit length for aligner.")
	flag.Float64Var(&dpMinId, "dpid", 0, "Minimum hit identity for aligner.")
	flag.IntVar(&tubeOffset, "tubeoffset", 0, "Tube offset - 0 indicate autotune.")

	flag.StringVar(&tmpDir, "tmp", "", "Path for temporary files.")
	flag.IntVar(&tmpChunk, "chunk", 1<<20, "Chunk size for morass.")
	flag.BoolVar(&tmpConcurrent, "tmpcon", false, "Process morass concurrently.")

	flag.IntVar(&threads, "threads", 1, "Number of threads to use for alignment.")
	flag.Uint64Var(&maxMem, "mem", 0, "Maximum nominal memory - 0 indicates unlimited.")

	flag.BoolVar(&logToFile, "log", false, "Log to file.")
	flag.BoolVar(&debug, "debug", false, "Include file names/lines in log.")
	flag.BoolVar(&verbose, "v", false, "Log additional information.")

	flag.StringVar(&cpuprofile, "cpuprofile", "", "write cpu profile to this file.")

	help := flag.Bool("help", false, "Print this help message.")

	flag.Parse()

	if *help {
		flag.Usage()
		os.Exit(1)
	}

	if maxMem == 0 {
		mem = nil
	} else {
		*mem = uintptr(maxMem)
	}

	runtime.GOMAXPROCS(threads)
}

func initLog(fileName string) {
	var w io.Writer = os.Stderr
	if fileName != "" {
		if file, err := os.Create(fileName); err == nil {
			fmt.Fprintln(file, strings.Join(os.Args, " "))
			w = io.MultiWriter(os.Stderr, file)
		} else {
			fmt.Printf("Error: Could not open log file: %v", err)
			os.Exit(1)
		}
	}

	logger = log.New(w, fmt.Sprintf("%s:", filepath.Base(os.Args[0])), log.Flags())
	if debug {
		logger.SetFlags(log.Flags() | log.Lshortfile)
	}
}

func main() {
	if cpuprofile != "" {
		profile, err := os.Create(cpuprofile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error: %v.", err)
			os.Exit(0)
		}
		fmt.Fprintf(os.Stderr, "Writing CPU profile data to %s\n", cpuprofile)
		pprof.StartCPUProfile(profile)
		defer pprof.StopCPUProfile()
	}

	if logToFile {
		initLog("krishna-" + time.Now().Format(timeFormat) + "-" + strconv.Itoa(pid) + ".log")
	} else {
		initLog("")
	}

	logger.Println(os.Args)
	var target, query *seq.Seq
	if targetName != "" {
		target = packSequence(targetName)
	} else {
		logger.Fatalln("No target provided.")
	}

	var writer *pals.Writer
	if outFile == "" {
		writer = pals.NewWriter(os.Stdout, 2, 60, false)
	} else {
		var err error
		writer, err = pals.NewWriterName(outFile, 2, 60, false)
		if err != nil {
			log.Fatalf("Could not open output file: %v", err)
		}
	}
	defer writer.Close()

	if !selfCompare {
		if queryName != "" {
			query = packSequence(queryName)
		} else {
			logger.Fatalln("No query provided in non-self comparison.")
		}
	} else {
		query = target
	}

	if maxK > 0 {
		pals.MaxKmerLen = maxK
	}

	m, err := morass.New(filter.FilterHit{}, "krishna_"+strconv.Itoa(pid)+"_", tmpDir, tmpChunk, tmpConcurrent)
	if err != nil {
		logger.Fatalf("Error: %v", err)
	}
	pa := pals.New(target, query, selfCompare, m, threads, tubeOffset, mem, logger)

	if err = pa.Optimise(minHitLen, minId); err != nil {
		logger.Fatalf("Error: %v", err)
	}
	if dpMinHitLen != 0 && dpMinId != 0 {
		pa.DPParams.MinHitLength = minHitLen
		pa.DPParams.MinId = minId
	}
	logger.Printf("Using filter parameters:")
	logger.Printf("\tWordSize = %d", pa.FilterParams.WordSize)
	logger.Printf("\tMinMatch = %d", pa.FilterParams.MinMatch)
	logger.Printf("\tMaxError = %d", pa.FilterParams.MaxError)
	logger.Printf("\tTubeOffset = %d", pa.FilterParams.TubeOffset)
	logger.Printf("\tAvg List Length = %.3f", pa.AvgIndexListLength(pa.FilterParams))
	logger.Printf("Using dynamic programming parameters:")
	logger.Printf("\tMinLen = %d", pa.DPParams.MinHitLength)
	logger.Printf("\tMinID = %.1f%%", pa.DPParams.MinId*100)
	logger.Printf("Estimated minimum memory required = %dMiB", pa.MemRequired(pa.FilterParams)/(1<<20))
	logger.Printf("Building index for %s", target.ID)

	if err = pa.BuildIndex(); err != nil {
		logger.Fatalf("Error: %v", err)
	}

	both := !sameStrand
	for _, comp := range [...]bool{false, true} {
		if comp {
			logger.Println("Working on complementary strands")
		} else {
			logger.Println("Working on self strand")
		}
		if both || !comp {
			hits, err := pa.Align(comp)
			if err != nil {
				logger.Fatalf("Error: %v", err)
			}

			logger.Println("Writing results")
			n, err := WriteDPHits(writer, target, query, hits, comp)
			if err != nil {
				logger.Fatalf("Error: %v.", err)
			}
			logger.Printf("Wrote hits (%v bytes)", n)
		}
	}

	pa.CleanUp()
}
