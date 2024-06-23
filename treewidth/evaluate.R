readDataFromFile <- function(fileName) {
  # Read all lines from the file
  lines <- readLines(fileName)
  
  # Convert each line to an integer
  int_values <- as.integer(lines)
  
  # Return the list of integer values
  return(int_values)
}

scatterplot_lists <- function(list1, list2, h1, h2, setting, filename) {
  # Check if the lengths of the lists are equal
  if (length(list1) != length(list2)) {
    stop("The lists must have the same length.")
  }
  
  png(filename=filename, width = 800, height = 600)
  
  max_value_list1 = max(list1)
  max_value_list2 = max(list2)
  max_value = max(max_value_list1, max_value_list2) + 10
  min_value_list1 = min(list1)
  min_value_list2 = min(list2)
  min_value = min(min_value_list1, min_value_list2) - 10
  listLabels <- 1:100
  
  # Create the scatter plot for the first list
  plot(
    x = listLabels, 
    y = list1, 
    main = "Heuristics", 
    xlab = "Graphs", 
    ylab = "Tw", 
    pch = 19,       # Solid circle for points
    col = "blue",   # Color of the points for List1
    xlim = range(c(1, 100)),
    ylim = range(c(min_value, max_value))
  )
  
  # Add points for the second list in a different color
  points(
    x = listLabels, 
    y = list2, 
    pch = 19,       # Solid circle for points
    col = "red"     # Color of the points for List2
  )
  
  # Add a legend
  legend("topright", legend = c(h1, h2, setting), col = c("blue", "red", "black"), pch = 19)

  dev.off()
}

plot_statistics <- function(list1, list2, list3, h1, h2, h3, setting, filename) {
  # Check if the lengths of the lists are non-zero
  if (length(list1) == 0 || length(list2) == 0 || length(list3) == 0) {
    stop("All lists must have non-zero length.")
  }
  
  # Calculate the average and median for each list
  avg1 <- mean(list1)
  med1 <- median(list1)
  avg2 <- mean(list2)
  med2 <- median(list2)
  avg3 <- mean(list3)
  med3 <- median(list3)
  
  # Combine the statistics into a data frame
  data <- data.frame(
    List = factor(rep(c("List1", "List2", "List3"), each = 2)),
    Statistic = rep(c("Average", "Median"), 3),
    Value = c(avg1, med1, avg2, med2, avg3, med3)
  )
  
  # Open a PNG device to save the plot
  png(filename, width = 800, height = 600)
  
  # Create the barplot
  barplot(
    height = data$Value, 
    beside = TRUE, 
    names.arg = data$Statistic, 
    col = c("blue","blue", "red","red", "green","green"), 
    args.legend = list(title = "Lists", x = "topright"),
    main = "Average and Median", 
    xlab = "Statistic", 
    ylab = "Value"
  )

  # Add a legend
  legend("topleft", legend = c(h1, h2, h3, setting), col = c("blue", "red", "green", "black"), pch = 19)  
  
  # Close the PNG device
  dev.off()
}


main <- function() {
  max_card_10_0.25 <- readDataFromFile("results/max_card_10_0.250000_results.dat")
  max_card_10_0.50 <- readDataFromFile("results/max_card_10_0.500000_results.dat")
  max_card_10_0.75 <- readDataFromFile("results/max_card_10_0.750000_results.dat")
  max_card_100_0.25 <- readDataFromFile("results/max_card_100_0.250000_results.dat")
  max_card_100_0.50 <- readDataFromFile("results/max_card_100_0.500000_results.dat")
  max_card_100_0.75 <- readDataFromFile("results/max_card_100_0.750000_results.dat")
  min_deg_10_0.25 <- readDataFromFile("results/min_deg_10_0.250000_results.dat")
  min_deg_10_0.50 <- readDataFromFile("results/min_deg_10_0.500000_results.dat")
  min_deg_10_0.75 <- readDataFromFile("results/min_deg_10_0.750000_results.dat")
  min_deg_100_0.25 <- readDataFromFile("results/min_deg_100_0.250000_results.dat")
  min_deg_100_0.50 <- readDataFromFile("results/min_deg_100_0.500000_results.dat")
  min_deg_100_0.75 <- readDataFromFile("results/min_deg_100_0.750000_results.dat")
  min_fill_in_10_0.25 <- readDataFromFile("results/min_fill_in_10_0.250000_results.dat")
  min_fill_in_10_0.50 <- readDataFromFile("results/min_fill_in_10_0.500000_results.dat")
  min_fill_in_10_0.75 <- readDataFromFile("results/min_fill_in_10_0.750000_results.dat")
  min_fill_in_100_0.25 <- readDataFromFile("results/min_fill_in_100_0.250000_results.dat")
  min_fill_in_100_0.50 <- readDataFromFile("results/min_fill_in_100_0.500000_results.dat")
  min_fill_in_100_0.75 <- readDataFromFile("results/min_fill_in_100_0.750000_results.dat")
  
  # plot the scatter plots
  # min deg vs min fill in
  scatterplot_lists(min_deg_10_0.25, min_fill_in_10_0.25, "min-degree", "min-fill-in", "n=10,p=0.25", "plots/min-deg_min-fill_n10_p25.png")
  scatterplot_lists(min_deg_10_0.50, min_fill_in_10_0.50, "min-degree", "min-fill-in", "n=10,p=0.50", "plots/min-deg_min-fill_n10_p50.png")
  scatterplot_lists(min_deg_10_0.75, min_fill_in_10_0.75, "min-degree", "min-fill-in", "n=10,p=0.75", "plots/min-deg_min-fill_n10_p75.png")
  scatterplot_lists(min_deg_100_0.25, min_fill_in_100_0.25, "min-degree", "min-fill-in", "n=100,p=0.25", "plots/min-deg_min-fill_n100_p25.png")
  scatterplot_lists(min_deg_100_0.50, min_fill_in_100_0.50, "min-degree", "min-fill-in", "n=100,p=0.50", "plots/min-deg_min-fill_n100_p50.png")
  scatterplot_lists(min_deg_100_0.75, min_fill_in_100_0.75, "min-degree", "min-fill-in", "n=100,p=0.75", "plots/min-deg_min-fill_n100_p75.png")
  
  # min def vs max card
  scatterplot_lists(min_deg_10_0.25, max_card_10_0.25, "min-degree", "max-card", "n=10,p=0.25", "plots/min-deg_max-card_n10_p25.png")
  scatterplot_lists(min_deg_10_0.50, max_card_10_0.50, "min-degree", "max-card", "n=10,p=0.50", "plots/min-deg_max-card_n10_p50.png")
  scatterplot_lists(min_deg_10_0.75, max_card_10_0.75, "min-degree", "max-card", "n=10,p=0.75", "plots/min-deg_max-card_n10_p75.png")
  scatterplot_lists(min_deg_100_0.25, max_card_100_0.25, "min-degree", "max-card", "n=100,p=0.25", "plots/min-deg_max-card_n100_p25.png")
  scatterplot_lists(min_deg_100_0.50, max_card_100_0.50, "min-degree", "max-card", "n=100,p=0.50", "plots/min-deg_max-card_n100_p50.png")
  scatterplot_lists(min_deg_100_0.75, max_card_100_0.75, "min-degree", "max-card", "n=100,p=0.75", "plots/min-deg_max-card_n100_p75.png")
  
  # min fill in vs max card
  scatterplot_lists(min_fill_in_10_0.25, max_card_10_0.25, "min-fill-in", "max-card", "n=10,p=0.25", "plots/min-fill_max-card_n10_p25.png")
  scatterplot_lists(min_fill_in_10_0.50, max_card_10_0.50, "min-fill-in", "max-card", "n=10,p=0.50", "plots/min-fill_max-card_n10_p50.png")
  scatterplot_lists(min_fill_in_10_0.75, max_card_10_0.75, "min-fill-in", "max-card", "n=10,p=0.75", "plots/min-fill_max-card_n10_p75.png")
  scatterplot_lists(min_fill_in_100_0.25, max_card_100_0.25, "min-fill-in", "max-card", "n=100,p=0.25", "plots/min-fill_max-card_n100_p25.png")
  scatterplot_lists(min_fill_in_100_0.50, max_card_100_0.50, "min-fill-in", "max-card", "n=100,p=0.50", "plots/min-fill_max-card_n100_p50.png")
  scatterplot_lists(min_fill_in_100_0.75, max_card_100_0.75, "min-fill-in", "max-card", "n=100,p=0.75", "plots/min-fill_max-card_n100_p75.png")


  # plot avg and median
  plot_statistics(min_deg_10_0.25, min_fill_in_10_0.25, max_card_10_0.25, "min-degree", "min-fill", "max-card", "n=10,p=0.25", "plots/avg_median_n10_p25.png")
  plot_statistics(min_deg_10_0.50, min_fill_in_10_0.50, max_card_10_0.50, "min-degree", "min-fill", "max-card", "n=10,p=0.50", "plots/avg_median_n10_p50.png")
  plot_statistics(min_deg_10_0.75, min_fill_in_10_0.75, max_card_10_0.75, "min-degree", "min-fill", "max-card", "n=10,p=0.75", "plots/avg_median_n10_p75.png")
  plot_statistics(min_deg_100_0.25, min_fill_in_100_0.25, max_card_100_0.25, "min-degree", "min-fill", "max-card", "n=100,p=0.25", "plots/avg_median_n100_p25.png")
  plot_statistics(min_deg_100_0.50, min_fill_in_100_0.50, max_card_100_0.50, "min-degree", "min-fill", "max-card", "n=100,p=0.50", "plots/avg_median_n100_p50.png")
  plot_statistics(min_deg_100_0.75, min_fill_in_100_0.75, max_card_100_0.75, "min-degree", "min-fill", "max-card", "n=100,p=0.75", "plots/avg_median_n100_p75.png")
 
  min_deg_total_results =  c(min_deg_10_0.25, min_deg_10_0.50, min_deg_10_0.75, min_deg_100_0.25,min_deg_100_0.50,min_deg_100_0.75)
  min_fill_in_total_results =  c(min_fill_in_10_0.25, min_fill_in_10_0.50, min_fill_in_10_0.75, min_fill_in_100_0.25,min_fill_in_100_0.50,min_fill_in_100_0.75)
  max_card_total_results =  c(max_card_10_0.25, max_card_10_0.50, max_card_10_0.75, max_card_100_0.25,max_card_100_0.50,max_card_100_0.75)
  plot_statistics(min_deg_total_results, min_fill_in_total_results, max_card_total_results, "min-degree", "min-fill", "max-card", "overall", "plots/avg_median_overall.png")
}

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

main()

