// Package to read and write FASTA format files
package fasta

// Copyright ©2011 Dan Kortschak <dan.kortschak@adelaide.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

import (
	"fmt"
	"bufio"
	"bytes"
	"github.com/kortschak/BioGo/bio"
	"github.com/kortschak/BioGo/seq"
	"github.com/kortschak/BioGo/util"
	"io"
	"os"
)

var (
	IDPrefix = []byte(">") // default delimiters
	SeqPrefix = []byte("") // default delimiters
)

// Fasta sequence format reader type.
type Reader struct {
	f         io.ReadCloser
	r         *bufio.Reader
	IDPrefix  []byte
	SeqPrefix []byte
	last      []byte
}

// Returns a new fasta format reader using f.
func NewReader(f io.ReadCloser) *Reader {
	return &Reader{
		f:         f,
		r:         bufio.NewReader(f),
		IDPrefix:  IDPrefix, 
		SeqPrefix: SeqPrefix,
		last:      nil,
	}
}

// Returns a new fasta format reader using a filename.
func NewReaderName(name string) (r *Reader, err error) {
	var f *os.File
	if f, err = os.Open(name); err != nil {
		return
	}
	return NewReader(f), nil
}

// Read a single sequence and return it or an error.
func (self *Reader) Read() (sequence *seq.Seq, err error) {
	var label, body []byte
	for {
		read, err := self.r.ReadBytes('>')
		if len(read) > 1 {
			// sanitize newlines
			read = bytes.Replace(read, []byte("\r\n"), []byte("\n"), -1)
			read = bytes.Replace(read, []byte("\n\r"), []byte("\n"), -1)
			read = bytes.Replace(read, []byte("\r"), []byte("\n"), -1)
					
			lines := bytes.Split(read, []byte("\n"))
			if len(lines) > 1 {
				label = lines[0]
				body = bytes.Join(lines[1:len(lines)-1], []byte(""))
			}
			break 
		} else if err != nil {
			return nil, io.EOF
		}
	}
	if len(label) > 0 && len(body) > 0 {
		sequence = seq.New(string(label), body, nil)
	} else {
		return nil, bio.NewError("Invalid fasta entry", 0, nil)
	}
	return
}

// Rewind the reader.
func (self *Reader) Rewind() (err error) {
	if s, ok := self.f.(io.Seeker); ok {
		self.last = nil
		_, err = s.Seek(0, 0)
		self.r = bufio.NewReader(self.f)
	} else {
		err = bio.NewError("Not a Seeker", 0, self)
	}
	return
}

// Close the reader.
func (self *Reader) Close() (err error) {
	return self.f.Close()
}

// Fasta sequence format writer type.
type Writer struct {
	f         io.WriteCloser
	w         *bufio.Writer
	IDPrefix  string
	SeqPrefix string
	Width     int
}

// Returns a new fasta format writer using f.
func NewWriter(f io.WriteCloser, width int) *Writer {
	return &Writer{
		f:         f,
		w:         bufio.NewWriter(f),
		IDPrefix:  ">", // default delimiters
		SeqPrefix: "",  // default delimiters
		Width:     width,
	}
}

// Returns a new fasta format writer using a filename, truncating any existing file.
// If appending is required use NewWriter and os.OpenFile.
func NewWriterName(name string, width int) (w *Writer, err error) {
	var f *os.File
	if f, err = os.Create(name); err != nil {
		return
	}
	return NewWriter(f, width), nil
}

// Write a single sequence and return the number of bytes written and any error.
func (self *Writer) Write(s *seq.Seq) (n int, err error) {
	n, err = self.w.Write([]byte(Format(s, self.Width)))
	return
}

// Format a single sequence into fasta string
func Format(s *seq.Seq, width int) (sequence string) {
	sequence = fmt.Sprintf("%s%s\n", IDPrefix, s.ID)
	if width > 0 {
		for i := 0; i*width <= s.Len(); i++ {
			endLinePos := util.Min(width*(i+1), s.Len())
			sequence = fmt.Sprintf("%s%s%s\n", sequence, SeqPrefix, string(s.Seq[width*i:endLinePos]))
		}
	} else {
		sequence = fmt.Sprintf("%s%s%s\n", sequence, SeqPrefix, s.Seq)
	}
	return 
}

// Flush the writer.
func (self *Writer) Flush() error {
	return self.w.Flush()
}

// Close the writer, flushing any unwritten sequence.
func (self *Writer) Close() (err error) {
	if err = self.w.Flush(); err != nil {
		return
	}
	return self.f.Close()
}
