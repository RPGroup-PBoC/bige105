# `docs`
This folder contains subfolders for each homework and the corresponding
solutions as `.tex` files and. The structure of each subfile should be as
follows.

```
docs
 |
 |----homework_1/
        |
        |---- figs/
        |      |
        |      |--- fig1.pdf
        |      |--- fig2.pdf
        |---- code/
        |      |
        |      |--- fig1.py
        |      |--- fig2.py
        |---- main.tex
        |---- prob1.tex
        |---- prob2.tex
```

With this commonly shared structure, we can each work on the homework files
without having to simultaneously edit a single file. It also makes
troubleshooting easier. 

## Contents 

The `main.tex` file in each homework should have the same preamble, but with
different titles and quotes, if desired. Within the `\begin{document} ...
\end{document}` environments, each problem will be included using the
`\input{prob1.tex}` or `\include{prob1.tex}` commands. 

Solutions can be written within each problem portion between `\begin{solution}
... \end{solution}` environments. A key variable in the `main.tex` file can be
toggled to make the solutions visible/invisible. 

## Working on problems
While working on problems or assigning the work to others, try to keep track of
this using the GitHub issue tracker. Each homework will have its own issue which
will be closed when the homework is completed, meaning solutions and the main
content. 

