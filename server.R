shinyServer(function(input, output){ 
  # source data & app tools from base
  for (file in list.files(c(file.path(r_path,"base")),
                          pattern="\\.(r|R)$", full.names = TRUE))
    source(file, encoding = r_encoding, local = TRUE)
  
})
