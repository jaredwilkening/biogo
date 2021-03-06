package tree

// Copyright ©2011 Dan Kortschak <dan.kortschak@adelaide.edu.au>
// Derived from PyCogent tree package Copyright ©2007-2011, The Cogent Project, under GPL2 or greater
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

const (
	PreOrder = iota
	PostOrder
	PrePostOrder
	LevelOrder
)

type Node struct {
	Name       string
	support    float32
	length     float32
	parameters map[string]float32
	children   NodeList
	parent     *Node
	tree       *Tree
}

func NewNode(name string, support, length float32) *Node {
	return &Node{
		Name:    name,
		support: support,
		length:  length,
	}
}

func (self *Node) AddNode(n *Node) {
	self.children.Push(n)
	n.parent = self
	n.tree = self.tree
	self.tree.nodes.Push(n)
	for _, c := range n.children.nodeList {
		self.tree.nodes.Push(c)
	}
}

func (self *Node) NodeIterator(order byte, includeSelf bool, reaper <-chan struct{}) (c chan *Node) {
	c = make(chan *Node)
	switch order {
	case PreOrder:
		self.PreOrder(includeSelf, c, reaper)
	case PostOrder:
		self.PostOrder(includeSelf, c, reaper)
	case PrePostOrder:
		self.PrePostOrder(includeSelf, c, reaper)
	case LevelOrder:
		self.LevelOrder(includeSelf, c, reaper)
	}
	return
}

func (self *Node) PreOrder(includeSelf bool, c chan *Node, reaper <-chan struct{}) {
	go func() {
		defer close(c)

		var (
			i    int
			this *Node
		)

		stack := []*Node{self}
		for len(stack) > 0 {
			i = len(stack) - 1
			this, stack = stack[i], stack[:i]
			if this == self || includeSelf {
				select {
				case c <- this:
				case <-reaper:
					return
				}
			}
			if this.children.Len() > 0 {
				stack = append(stack, append(this.children.nodeList, stack...)...)
			}
		}
	}()
}

func (self *Node) PostOrder(includeSelf bool, c chan *Node, reaper <-chan struct{}) {
	go func() {
		defer close(c)

		var (
			index       int
			this, child *Node
		)

		childIndex := []int{0}
		this = self
		for {
			index = childIndex[len(childIndex)-1]
			if index < len(this.children.nodeList) {
				child = this.children.nodeList[index]
				if len(child.children.nodeList) > 0 {
					childIndex = append(childIndex, 0)
					this = child
					index = 0
				} else {
					select {
					case c <- child:
					case <-reaper:
						return
					}
					childIndex[len(childIndex)-1]++
				}
			} else {
				if includeSelf || this != self {
					select {
					case c <- this:
					case <-reaper:
						return
					}
				}
				if this == self {
					break
				}
				this = this.parent
				childIndex = childIndex[:len(childIndex)-1]
				childIndex[len(childIndex)-1]++
			}
		}
	}()
}

func (self *Node) PrePostOrder(includeSelf bool, c chan *Node, reaper <-chan struct{}) {
	go func() {
		defer close(c)

		if self.children.Len() < 1 {
			if includeSelf {
				select {
				case c <- self:
				case <-reaper:
					return
				}
			}
		} else {
			var (
				this, child *Node
				i, index    int
			)

			childIndex := []int{0}
			this = self
			for {
				index = childIndex[len(childIndex)-1]
				if index < 1 {
					if this != self || includeSelf {
						select {
						case c <- self:
						case <-reaper:
							return
						}
					}
				}
				if index < len(this.children.nodeList) {
					child = this.children.nodeList[index]
					if len(child.children.nodeList) > 0 {
						childIndex = append(childIndex, 0)
						this = child
						index = 0
					} else {
						select {
						case c <- child:
						case <-reaper:
							return
						}
						childIndex[len(childIndex)-1]++
					}
				} else {
					if includeSelf || this != self {
						select {
						case c <- this:
						case <-reaper:
							return
						}
					}
					if this == self {
						break
					}
					this = this.parent
					i = len(childIndex) - 1
					childIndex = childIndex[:i]
					childIndex[i]++
				}
			}
		}
	}()
}

func (self *Node) LevelOrder(includeSelf bool, c chan *Node, reaper <-chan struct{}) {
	go func() {
		defer close(c)

		var (
			this, child *Node
		)

		queue := []*Node{self}

		for len(queue) > 0 {
			this = queue[0]
			queue = append(queue[:0], queue[1:]...)
			if this != self || includeSelf {
				select {
				case c <- this:
				case <-reaper:
					return
				}
			}
			if len(this.children.nodeList) > 0 {
				for _, child = range this.children.nodeList {
					queue = append(queue, child)
				}
			}
		}
	}()
}

func (self *Node) Nodes(order byte, includeSelf bool) (nodes NodeList) {
	for n := range self.NodeIterator(order, includeSelf, nil) {
		nodes.Push(n)
	}
	return
}

func (self *Node) InternalNodeIterator(includeSelf bool, reaper <-chan struct{}) (c chan *Node) {
	c = make(chan *Node)
	go func() {
		defer close(c)

		r := make(chan struct{})
		for n := range self.NodeIterator(PreOrder, includeSelf, r) {
			if n.children.Len() > 0 {
				select {
				case c <- n:
				case <-reaper:
					close(r)
					return
				}
			}
		}
	}()
	return
}

func (self *Node) InternalNodes(includeSelf bool) (internals NodeList) {
	for n := range self.LeafIterator(includeSelf, nil) {
		internals.Push(n)
	}
	return
}

func (self *Node) LeafIterator(includeSelf bool, reaper <-chan struct{}) (c chan *Node) {
	c = make(chan *Node)
	go func() {
		defer close(c)

		if self.children.Len() < 1 {
			if includeSelf {
				select {
				case c <- self:
				case <-reaper:
					return
				}
			}
		} else {
			var (
				i    int
				this *Node
			)

			stack := []*Node{self}
			for len(stack) > 0 {
				i = len(stack) - 1
				this, stack = stack[i], stack[:i]
				if len(this.children.nodeList) > 0 {
					stack = append(stack, append(this.children.nodeList, stack...)...)
				} else {
					select {
					case c <- this:
					case <-reaper:
						return
					}
				}
			}
		}
	}()
	return
}

func (self *Node) Leaves(includeSelf bool) (leaves NodeList) {
	for n := range self.LeafIterator(includeSelf, nil) {
		leaves.Push(n)
	}
	return
}
