library(bio3d)

# TODO: seq_len emits a warning about length.out:
# "Warning: first element used of 'length.out' argument"

map_scores_pdb <- function(input_pdb, mapping_scores, field, selection = NULL) {
  
  if (is.null(selection)) {
    
    selection = atom.select(input_pdb, "protein")
  }
  
  output_pdb = trim.pdb(input_pdb, selection)
  
  for (i in seq_len(dim(output_pdb$atom[1]))) {
    
    if (output_pdb$atom[i,]$resno > 0) {
      
      n = as.character(output_pdb$atom[i,]$resno)
      j = which(mapping_scores['pos'] == n)

      if (length(j) == 0) {
        score = 0
        
      } else {
        score = mapping_scores[j, field][[1]]
      }
      
      if (!is.na(score)) {
        
        output_pdb$atom[i,]$b = score
        
      } else {
        
        output_pdb$atom[i,]$b = 0
        
      }
    } else {
      output_pdb$atom[i,]$b = 0
    }
  }
  
  return(output_pdb)
}