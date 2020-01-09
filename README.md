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

## Deployment

This was deployed on an AWS EC instance running Ubuntu 16.04 LTS or 18.04 LTS. To install R, Shiny, and Shiny Server, I used this guide: https://rstudio.com/products/shiny/download-server/ubuntu/

After installing R, make sure to install the required packages (ggplot2, shinyBS, gmp, arrangements) using e.g.
`sudo su - -c "R -e \"install.packages('arrangements', repos='https://cran.rstudio.com/')\""`.

Make sure to configure Shiny Server to listen on port 80 in `/etc/shiny-server/shiny-server.conf`. Example configuration:
```
# Instruct Shiny Server to run applications as the user "shiny"
run_as shiny;

# Define a server that listens on port 3838
server {
  listen 80;

  server_name n1sim.hd2i.org;

  # Define a location at the base URL
  location / {

    # Host the directory of Shiny Apps stored in this directory
    # app_dir /srv/shiny-server-n1sim/shiny;
    app_dir /home/mikejones/n1-simulator/shiny;

    # Log all Shiny output to files in this directory
    log_dir /var/log/shiny-server-n1sim/;

    # When a user visits the base URL rather than a particular application,
    # an index of the applications available in this directory will be shown.
    # directory_index on;
  }
}
```

