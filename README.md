# n1-simulator

Simulator for exploring N-of-1 studies.
Accompanies _Percha, B., Baskerville, E. B., Johnson, M., Dudley, J. T., & Zimmerman, N. (2019). [Designing Robust N-of-1 Studies for Precision Medicine: Simulation Study and Design Recommendations](https://www.jmir.org/2019/4/e12641/). Journal of Medical Internet Research, 21(4), e12641._ which contains a full description of the methods.


## Installation

Only one non-standard package is required, `nlme`:

```{r}
install.packages('nlme')
```

All functions are in a single R file, which you can copy into your project and source:

```{r}
source('n1-simulator.R')
```

The examples use `ggplot2` for plotting:

```{r}
install.packages('ggplot2')
```

## Usage

See `examples/` in this repository.

## Shiny web application

To try out the simulator using a web user interface, visit http://n1sim.hd2i.org. See `shiny/` in this repository for the code.
