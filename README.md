# prodigy: A library for manipulation distributions through generating functions

The main features are
- python bindings
- a purely string-based interface, no need to mess around with data structures for algebraic expression and the like
- support for parameters (can also be used as meta symbols in second-order PGFs)

## Prerequisites

GiNaC must be installed on your system.
On macOS with homebrew type
```brew install ginac```

## Installation

Activate the virtual environment where prodigy should be installed, then type

```pip install .```

## Usage

```
import prodigy
dist = prodigy.geometric("x", "1/2")
e = dist.E("x")  # e becomes 2
'''
