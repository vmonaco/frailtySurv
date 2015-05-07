<!-- README.md is generated from README.Rmd. Please edit that file -->
frailtyr
--------

frailtyr is a library for fitting shared frailty models with support for general frailty distributions. **This is currently a work in progress**, part of Google Summer of Code 2015.

Installation
------------

To get started, you can install the package from github using `devtools`.

``` r
library(devtools)
install_github("vmonaco/frailtyr")
```

Clone and build
---------------

Instead of installing from github, you can clone and build the project from source. To create a project in RStudio from this repository, you must have both RStudio and git installed.

-   In RStudio, go to File -\> New Project -\> Version Control -\> Git
-   Name the project (eg. frailtyr), choose location, and specify the Repository URL as `https://github.com/vmonaco/frailtyr`

More information about [RStudio and git](https://support.rstudio.com/hc/en-us/articles/200532077-Version-Control-with-Git-and-SVN)

Also see the chapter in the [R packages book](http://r-pkgs.had.co.nz/git.html)

### Build everything

To load the package (simulate installing and loading with `library(frailty)`):

``` r
devtools::load_all()
```

or shortuct: **Ctrl+Shift+L**

Use this command to build a source package (binary = FALSE by default)

``` r
devtools::build()
```

Other usefull commands and shortcuts:

-   Build and reload everything: **Ctrl+Shift+B**
-   Check the package, `devtools::check()`: **Ctrl+Shift+E**
-   Run all tests, `devtools::test()`: **Ctrl+Shift+T**

### Build vignettes

Build the vignettes with the command:

``` r
devtools::build_vignettes()
```

Compiled vignettes will reside in `inst/doc`. In particular, see frailtyr.pdf

### Run tests

frailtyr uses the `testthat` library for tests. Test files are all located in `tests/testthat`. To run all tests, use the command `devtools::test()` or shortcut **Ctrl+Shift+T**

### git pull

Use RStudio to pull the latest code if you've already cloned the repository. In the "git" pane, look for the Pull button. You can also view the commit history from there and see the diff between modified files, if any.
