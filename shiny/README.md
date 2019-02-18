# n1-simulator Shiny app
This is a web application that adds a user interface and dynamic data visualization to the [simulator code](https://github.com/HD2i/n1-simulator) in this repository.

## Use
The web application can be used directly from http://n1sim.hd2i.org.

## Local installation
To run a local instance of this R Shiny app:
1. Clone this repository.
2. Install [R Studio Desktop](https://www.rstudio.com/products/rstudio/download/#download).
3. Open R Studio Desktop and run the following commands from the R Studio Desktop console to install required packages:
   * install.packages("shiny")
   * install.packages("ggplot2")
   * install.packages("arrangements")
   * install.packages("shinyBS")
4. In R Studio Desktop, open the file n1-simulator/shiny/app.R.
5. Click the button "Run App" to run the Shiny application.
6. Open the localhost link that is displayed in the R Studio Desktop console using a web browser (e.g. http://127.0.0.1:6025).

