library("sevenbridges")
a <- Auth(platform = "cavatica", token = "d8554e151de14d9b998af569994f4bfb")
p <- a$project(id = "gaonkark/sv-test")

files <- p$file (
  name = "*.maf",
  metadata = list(
    sample_type = "Normal"
    ), complete=TRUE)

filenames <- lapply(files, function(x) x$name)

generate_CPV_maf_task<-function(x){
  input_file<-x[1]
  tsk_name=paste("CPV_maf",input_file,date(),sep="-")
  tsk <- p$task_add(
    name = tsk_name,
    app = "gaonkark/sv-test/maf-filter",
    inputs = list(
      input_maf=input_file,
      input_genes="predisposition-genes.txt",
      threshold = 0.001
    )
  )
  tsk$run()
}

lapply(filenames, function(x) generate_CPV_maf_task(x))



