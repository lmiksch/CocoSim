# Cocopaths - A Compiler for Cotranscriptional Folding Pathways

[![codecov](https://codecov.io/gh/lmiksch/CocoSim/graph/badge.svg?token=ROMQLY0I9T)](https://codecov.io/gh/lmiksch/CocoSim)

Cocopaths is a compiler designed to analyze cotranscriptional folding pathways. This repository includes **Cocosim**, a tool to simulate folding processes interactively or via input files.

## Getting Started

Clone the repository and install the package:

    # Clone this repository
    git clone https://github.com/lmiksch/cocosim

    # Change your directory to the cloned repository
    cd cocosim

## Installation

Install the package using pip:

- **Standard Installation:**

      pip install .

- **Development Installation (editable mode):**

      pip install -e .

## Usage

After installation, you can use **Cocosim** in different ways.

### Running Tests

To ensure all tests pass, run:

      pytest

### Interactive Mode

Running **Cocosim** without any options will prompt you to enter a domain-level sequence:

      cocosim

You'll see a prompt like:

      Please enter a domain level sequence:

If you do not have an aCFP (a specific folding path) available, you can skip this step by entering `$`.

### Input File Mode

You can also provide input in `.pil` format. An example input file is available in the `examples` folder. To run a simulation using the example file, use:

      cocosim -i examples/example_input.pil

This command will process the input file and simulate the folding pathway based on its content.

### Using the -l Flag

The `-l` flag runs Cocosim in a mode that displays logical domain pairing information along with the folding simulation. When using `-l`, Cocosim prompts you for a domain-level sequence and then interactively collects your folding path inputs in dot-bracket notation.

To run Cocosim with the `-l` flag, execute:

      cocosim -l

A sample session might look like this:

    (test-cocosim) miksch@romulus:~/Documents/CocoSim$ cocosim -l 
    Please enter a domain level sequence:
    a L0 b L0* a c a* L0 
    Given Domain Sequence: a L0 b L0* a c a* L0 

    Current Input: []
    Please input a folding path in dot-bracket annotation or use '$' to exit input and continue use 'r' to reset input:
    .

    Current Input: ['.']
    Please input a folding path in dot-bracket annotation or use '$' to exit input and continue use 'r' to reset input:
    ()

    Current Input: ['.', '()']
    Please input a folding path in dot-bracket annotation or use '$' to exit input and continue use 'r' to reset input:
    .()

    Current Input: ['.', '()', '.()']
    Please input a folding path in dot-bracket annotation or use '$' to exit input and continue use 'r' to reset input:
    $

    Final Input:
    ['.', '()', '.()']

    Cocosim
    Domain Level seq: a L0 b L0* a c a* L0 

    Transcription Step | Occupancy  |  Logic domain pairing | Structure     
      2   |   1.0000   |   .               |   a L0
      3   |   1.0000   |   .               |   a L0 b
      4   |   1.0000   |   ..              |   a L0 b L0*
      5   |   1.0000   |   ()              |   a L0( b ) a
      6   |   1.0000   |   ()              |   a L0( b ) a c
      7   |   1.0000   |   ()              |   a L0( b ) a c a*
      8   |   0.1095   |   ().            |   a L0( b ) a c a* L0
      8   |   0.8420   |   ().            |   a L0( b ) a( c ) L0
      8   |   0.0485   |   ().            |   a( L0( b ) a c ) L0

    END  |   1.0000   |   .()             |   a L0 b L0*( a( c ) )

In this mode, after you input your folding path (using dot-bracket notation), Cocosim displays the transcription steps along with:
- **Occupancy:** The extent of structure formation.
- **Logic domain pairing:** A representation of the pairing logic based on your input.
- **Structure:** The evolving domain-level structure.


## Dependencies

**Cocosim** depends on the following library:

- [ViennaRNA](https://www.tbi.univie.ac.at/RNA/)

Ensure that this dependency is installed and accessible in your environment.

## License

This project is licensed under the terms specified in the [LICENSE](LICENSE) file.
