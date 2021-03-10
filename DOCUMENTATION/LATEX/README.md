# The Martian Report 
**TMR** is a LaTeX template for written reports, written with sleek modern design in mind.

# Preview
| Cover Page | Section Heading | Body |
|:---:|:---:|:---:|
| [![Cover Page](/examples/cover-page.png?raw=true)](Example.pdf) | [![Section Heading](/examples/section-heading.png?raw=true)](Example.pdf) | [![Body](/examples/body.png?raw=true)](Example.pdf) |

# Usage
This package uses [**The MIT License**](https://opensource.org/licenses/MIT). 

It requires the [geometry](http://mirror.aarnet.edu.au/pub/CTAN/macros/latex/contrib/geometry/geometry.pdf), [fancyhdr](http://mirror.aarnet.edu.au/pub/CTAN/macros/latex/contrib/fancyhdr/fancyhdr.pdf), [parskip](http://mirror.aarnet.edu.au/pub/CTAN/macros/latex/contrib/parskip/parskip-doc.pdf), [xcolor](http://mirror.aarnet.edu.au/pub/CTAN/macros/latex/contrib/xcolor/xcolor.pdf), [hyperref](http://ftp.math.purdue.edu/mirrors/ctan.org/macros/latex/contrib/hyperref/doc/manual.pdf), [wrapfig](http://ftp.math.purdue.edu/mirrors/ctan.org/macros/latex/contrib/wrapfig/wrapfig-doc.pdf), [graphicx](http://mirrors.rit.edu/CTAN/macros/latex/required/graphics/grfguide.pdf), [tocloft](http://ctan.math.utah.edu/ctan/tex-archive/macros/latex/contrib/tocloft/tocloft.pdf), [xparse](http://ctan.math.utah.edu/ctan/tex-archive/macros/latex/contrib/l3packages/xparse.pdf) & [natbib](http://mirrors.ibiblio.org/CTAN/macros/latex/contrib/natbib/natbib.pdf) packages, all of which use [**The LATEX Project Public License**](http://www.latex-project.org/lppl) (LPPL).

# How to Use

## Values

| Value | Description | Default |
| --- | --- | --- |
| **\title**{*text*} | Document title string for cover page & footer) | "Concept Document" |
| **\author**{*text*} | Document author string (for cover page & footer) | blank |
| **\project**{*text*} | Document project name (optional) | blank |
| **\imagefolder**{*relative-path/*} | Folder to look for images in | . (here) |
| **\revision**{*value*} | Document revision text/number | "Revision 1.0" |
| **\footertext**{*left*}{*right*} | Text to show in footer | \@title & \@author |
| **\imageref**{*text*} | Text to show at end of document, used for referencing header images | blank |
| **\bibliographystyle**{*stylename*} | Style to format References in | abbrvnat |
| **\logo**{*image-filename.jpg*} | Image file to display on cover page | none (hidden) |

## Commands

|Command | Use |
| --- | --- |
| **\printhelp** | Adds page to PDF that shows all available TMR commands and help |
| **\drafting** | Hides \temp and \crit sections in text for final printing |
| **\nofootertext** | Equivalent to \footertext{}{} |
| **\sectionnumbers** | Displays numbers in section/subsection headings |
| **\subsectionnumbers** | Displays numbers in subsection/subsubsection headings |
| **\section**{*name*}{*image-filename.jpg*} | Same as \section but second value is image name to display above heading |
| **\temp**{*text*}	| Text to show as blue italics, hidden if \drafting is set |
| **\crit**{*text*}	| Text to show as red italics, hidden if \drafting is set |
| **\link**{*url*}{*optional text*} | Text to display as link to url, else just url |
| **\image**{*width*}{*image-filename.jpg*}{*optional text*} | Image to include as width=*width*\textwidth, with optional figure caption |
| **\begin{nicelist}**{*text*}... | Same as itemize but with *text* as title |
| **\begin{niceenum}**{*text*}... | Same as enumerate but with *text* as title |
| **\tipbox**{*text*}{*text*} | Makes box with first parameter as title and second as body text (use \tipbox*{}{} to wrap text) |