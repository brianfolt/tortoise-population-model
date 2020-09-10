

# Function 1: create known latent states z
known.state.ms <- function(ms, notseen){  
  # Requires two variables: cap history & code for not seen
  # notseen: label for 'not seen'
  state <- ms
  state[state==notseen] <- NA    # if not seen, put NA
  for (i in 1:dim(ms)[1]){       # for all individuals...
    m <- min(which
             
             (!is.na(state[i,])))  # m is minimum/first occasion that's not equal to NA for that indiv
    state[i,m] <- NA            # first capture occs'n is assigned NA
  }
  return(state)
}

# Function 2: create initial values and NA for first encounter
ms.init.z  <- function(ms, notseen){
  state <- ms					# Capture history named state
  state[state==notseen] <- NA	# If not seen, then NA
  
  for(i in 1:dim(state)[1]){	# For every individual (row)
    # If any observation is not NA through the 2nd to last occasion
    if(any(!is.na(state[i,1:dim(state)[2]-1]))){   
      # Apply max state value for that individual to final occasions
      state[i,dim(state)[2]] <- max(state[i,],na.rm=TRUE) 
      # Populate last column b/c as.interpolation needs at least two values
      state[i,] <- ceiling(na_interpolation(state[i,])) # Interpolate
    } # End-if
    m <- min(which(!is.na(ms[i,])))   # Identify the first occasion not NA
    state[i,1:m] <- NA	# Replace before and on first occasion with NA  
  } # End-if
  
  for(i in 1:dim(state)[1]){  # For each row in the capture history  
    #Identify situations where juveniles (1) were observed to transitioned 
    # to become adult males (3)
    # AMT - 10 sept 19 - changed this from 
    # if(state[i,dim(state)[2]-1] == 3 ...)
    # by looking at second to last occasion was missing some instances where ind was not seen as an adult
    # unti the last year
    if(state[i,dim(state)[2]]==3 & 2 %in% state[i,1:dim(state)[2]]){
      # If max state = 3 (male) and contains a previous state = 1 (juvenile)
      state[i,] = replace(state[i,], state[i,]==2, 3)
      # Force any instances of female interpolation (2) to be male (3) 
    } # End-if
  } # End-i-loop
  
  return(state)
} #function
