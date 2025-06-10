Thoughts:

Would this be a general and easy to use setup:
The main data structure that the user supplies is a 4-column data frame. 
First 2 columns are full names of the strains, other 2 columns are sequences for each strain.
We could allow other columns as long as the first 4 are the ones we need. 
Or we require specific names of these columns.

Main functionality is something like this

dist <- calculate_pepitope_dist(dataframe, append = TRUE/FALSE)

which returns either a single vector with distances, or if append is true, 
returns a N+1 column data frame with the original N ones and the distance measure at the end.

Other distance measures follow the same structure. 

For calculate_temporal_dist(), the function would parse the dates in the file names.
If those are not provided, it will give an error. 
Alternatively, the user can specify the name of 2 columns which contain the dates.






Can we structure code such that each R file contains a single function of the same name. 


