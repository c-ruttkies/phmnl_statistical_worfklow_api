#!/usr/bin/env Rscript
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)

# Taking the command line arguments
args <- commandArgs(trailingOnly = TRUE)

if(length(args)==0) {
  help.msg <- paste("", "Usage: ./run_statistical_workflow.r isatab2w4m=... study_path=... study_name=... ga_file_template=... output_path=... log_file=... galaxy_url=... galaxy_key=...",
    "",
    "",
    "Arguments:",
    "",
    "\tisatab2w4m - path to isatab2w4m script (https://github.com/workflow4metabolomics/mtbls-dwnld/tree/develop)",
    "",
    "\tstudy_path - local path containing MTBLS studies",
    "",
    "\tstudy_name - name of the study to be processed (e.g. MTBLS404)",
    "",
    "\tga_file_template - Galaxy workflow file for the statistical workflow",
    "",
    "\toutput_path - local path used to store the result files",
    "",
    "\tlog_file - path to csv file used for log messages",
    "",
    "\tgalaxy_url - url to Galaxy server on which to run the workflow",
    "",
    "\tgalaxy_key - API key to access the Galaxy server",
    "",
    "",
    "Description:",
    "",
    paste("\tRuns the statistical worfklow (https://github.com/phnmnl/phenomenal-h2020/wiki/Sacurine-statistical-workflow) using wft4galaxy.", 
      "\tMake sure that wft4galaxy-docker is available on your machine. Visit http://wft4galaxy.readthedocs.io/installation.html#id2 for installation instructions.",      "\tThe R package Risa needs to be installed. Visit https://bioconductor.org/packages/release/bioc/html/Risa.html for installation instructions.",
      sep="\n"),
    "",
    "",
    "Example: ",
    "",
    paste("\t./run_statistical_workflow.r isatab2w4m=\"mtbls-dwnld/isatab2w4m\"", 
      "\t\tstudy_path=\"studies\"",
      "\t\tstudy_name=\"MTBLS404\"",
      "\t\tga_file_template=\"template/w4m-sacurine-statistics.ga\"",
      "\t\toutput_path=\"output\"",
      "\t\tlog_file=\"MTBLS404_logs.csv\"",
      "\t\tgalaxy_url=\"http://193.62.54.91:30700\"",
      "\t\tgalaxy_key=\"ce30966564d6a42b22f951ca023081ed\"", sep=" \\\n"),
    "",
    sep = "\n"
  )
  message(help.msg)
  quit(save = "no", status = 1, runLast = TRUE)
}

sacurine.script <- NA
study.path <- NA
study.name <- NA
path.to.ga.template <- NA
output <- NA
log.file <- NA
galaxy.url <- NA
galaxy.key <- NA

for(arg in args)
{
  argCase<-strsplit(x = arg,split = "=")[[1]][1]
  value<-strsplit(x = arg,split = "=")[[1]][2]
  if(argCase=="isatab2w4m")
  {
    sacurine.script=as.character(value)
  }
  if(argCase=="study_path")
  {
    study.path=as.character(value)
  }
  if(argCase=="study_name")
  {
    study.name=as.character(value)
  }
  if(argCase=="ga_file_template")
  {
    path.to.ga.template=as.character(value)
  }
  if(argCase=="output_path")
  {
    output=as.character(value)
  }
  if(argCase=="log_file")
  {
    log.file=as.character(value)
  }
  if(argCase=="galaxy_url")
  {
    galaxy.url=as.character(value)
  }
  if(argCase=="galaxy_key")
  {
    galaxy.key=as.character(value)
  }
  
}

if(is.na(sacurine.script))stop("isatab2w4m needs to be specified!\n")
if(is.na(study.path))stop("study_path needs to be specified!\n")
if(is.na(study.name))stop("study_name needs to be specified!\n")
if(is.na(path.to.ga.template))stop("ga_file_template needs to be specified!\n")
if(is.na(output))stop("output_path needs to be specified!\n")
if(is.na(log.file))stop("log_file needs to be specified!\n")
if(is.na(galaxy.url))stop("galaxy_url needs to be specified!\n")
if(is.na(galaxy.key))stop("galaxy_key needs to be specified!\n")

if(!file.exists(sacurine.script))stop(paste("File path for isatab2w4m (", sacurine.script, ")", "not found.", sep=" "))
if(!dir.exists(study.path))stop(paste("Directory path for study_path (", study.path, ")", "not found.", sep=" "))
if(!dir.exists(paste(study.path, study.name, sep="/")))stop(paste("Study", study.name, "in", study.path, "not found. Check value of study_name.", sep=" "))
if(!dir.exists(output))stop(paste("Directory path for output_path (", output, ")", "not found.", sep=" "))
if(!file.exists(path.to.ga.template))stop(paste("File path for ga_file_template (", path.to.ga.template, ")", "not found.", sep=" "))
if(Sys.which("docker")=="")stop("docker seems not to be installed.")
if(system("wft4galaxy-docker --help", ignore.stdout = T, ignore.stderr = T) != 0)stop("The docker image for wft4galaxy-docker is not installed. Visit http://wft4galaxy.readthedocs.io/installation.html#id2 for instructions.")

library(Risa)

delete.dirs <- F

# do not modify this unless you know what you are doing
wft4galaxy.template.yaml <- paste("enable_logger: False
output_folder: \"RESULTS\"
# workflow tests
workflows:
 sacurine:
  file: \"w4m-sacurine-statistics.ga\"
  params:
   3:
    \"facC\": \"FACTOR_VALUE\"
   4:
    \"orthoI\": \"1\"
    \"predI\": \"1\"
    \"respC\": \"FACTOR_VALUE\"
    \"testL\": \"FALSE\"
   5:
    \"respC\": \"FACTOR_VALUE\"
  inputs:
   \"DataMatrix\": \"MATRIXDATA\"
   \"SampleMetadata\": \"SAMPLEDATA\"
   \"VariableMetadata\": \"VARIABLEDATA\"
  expected:
   Univariate_variableMetadata: \"expected/Univariate_variableMetadata.tsv\"
   Multivariate_sampleMetadata: \"expected/Multivariate_sampleMetadata.tsv\"
   Multivariate_variableMetadata: \"expected/Multivariate_variableMetadata.tsv\"
   Biosigner_variableMetadata: \"expected/Biosigner_variableMetadata.tsv\"", sep="")    

# retrieve modified factor name
get.sacurine.factor.name <- function(sample.file, variable.name) {
    data<-read.table(sample.file, header=T)
    regex.column.name<-paste(gsub("$",".*",gsub("\\s+", ".*", variable.name)), "$", sep="")
    found.indexes <- grep(regex.column.name, colnames(data))
    return(colnames(data)[found.indexes])
}
 
# checks whether assay file is MS type
check.assay.file <- function(assay.file) {
    data <- read.csv(assay.file, sep="\t")
    return(any(colnames(data)=="MS.Assay.Name"))
}
    
# function to check whether variable in sample file contains only 2 different values
check.sample.file <- function(sample.file, factor.name) {
    data <- read.table(sample.file, header=T)
    if(dim(data)[1] == 0) {
        return(paste("Error 6: Sample file data is empty: ", sample.file ,sep=""))
    }
    a <- c(length(unique(as.character(data[,factor.name]))))
    data[,factor.name]<-as.character(paste("c_",data[,factor.name],sep=""))
    write.table(data, file=sample.file, col.names=T, row.names=F,sep="\t")
    if(length(a) == 0 || a != 2) {
        return(paste("Error 7: Found class number unequal 2 for ", factor.name, " (",length(a), ")" ,sep=""))
    }
    return("")
}

# function to check whether variable in matrix file contains NA's
check.matrix.file <- function(matrix.file) {
    data <- read.table(matrix.file, header=T, stringsAsFactors=F)
    rownames(data) <- data[,1]
    data <- data[,c(-1)]
    a <- sapply(1:dim(data)[2], function(x) {
        if(!any(is.na(data[,x]))) {
        
        }
    })
    write.table(data, file=matrix.file, col.names=T, row.names=T, sep="\t")
    return("")
}
# generate wft4galaxy files
prepare.wft4galaxy.files <- function(factor.name, factor.folder, wft4galaxy.template.yaml, path.to.ga.template) {
    # create neccessary folders
    # folder for string results
    dir.create(paste(factor.folder,"results",sep="/"))
    # folder for storing expected results
    dir.create(paste(factor.folder,"expected",sep="/"))
    # prepare and write yaml file
    wft4galaxy.template.yaml.tmp <- gsub("RESULTS", "results", wft4galaxy.template.yaml)
    wft4galaxy.template.yaml.tmp <- gsub("FACTOR_VALUE", factor.name, wft4galaxy.template.yaml.tmp)
    wft4galaxy.template.yaml.tmp <- gsub("MATRIXDATA", "matrix-output.tsv", wft4galaxy.template.yaml.tmp)
    wft4galaxy.template.yaml.tmp <- gsub("SAMPLEDATA", "sample-output.tsv", wft4galaxy.template.yaml.tmp)
    wft4galaxy.template.yaml.tmp <- gsub("VARIABLEDATA", "variable-output.tsv", wft4galaxy.template.yaml.tmp)
    # write yaml file
    write(wft4galaxy.template.yaml.tmp, file=paste(factor.folder, "workflow.yaml", sep="/"))
    # create expected results
    file.create(paste(factor.folder,"expected/Univariate_variableMetadata.tsv",sep="/"))
    file.create(paste(factor.folder,"expected/Multivariate_sampleMetadata.tsv",sep="/"))
    file.create(paste(factor.folder,"expected/Multivariate_variableMetadata.tsv",sep="/"))
    file.create(paste(factor.folder,"expected/Biosigner_variableMetadata.tsv",sep="/"))
    # copy 
    file.copy(path.to.ga.template, paste(factor.folder, "w4m-sacurine-statistics.ga",sep="/"))
}

# run wft4galaxy in prepared folder
run.wft4galaxy <- function(factor.folder, galaxy.key, galaxy.url) {
    current.folder <- getwd()
    setwd(factor.folder)
    command <- paste(
        "export GALAXY_URL=",galaxy.url,"; ",
        "export GALAXY_API_KEY=",galaxy.key,"; ",
        "wft4galaxy-docker --debug runtest -f workflow.yaml --disable-cleanup --enable-logger",
        sep = ""
    )
    system(command)
    setwd(current.folder)
}

# check whether wft4galaxy run was successful
# otherwise return the error message
validate.wft4galaxy.run <- function(factor.folder, galaxy.url) {
    # get results folder
    result.folder.name <- paste(factor.folder, "results", sep="/")
    # search for log file
    log.file.index <- grep("\\.log$", dir(result.folder.name))
    if(length(log.file.index) != 0) {
        # check for error message
        lines <- readLines(paste(result.folder.name,dir(result.folder.name)[log.file.index],sep="/"))
        error.indexes <- grep("Runtime error:", lines)
        if(length(error.indexes) != 0) {
            # get first run id
            run.id <- gsub("\\s.*", "", gsub(".*Runtime error:\\s*", "", lines[error.indexes[1]]))
            download.result <- paste(factor.folder, "results", run.id, sep="/")
            url <- paste(galaxy.url,"/dataset/errors?id=", run.id, sep="")
            download.file(url, download.result, quiet = T)
            if(file.exists(download.result)) {
                lines <- readLines(download.result)
                # search for error messages
                error.started <- FALSE
                error.message <- ""
                for(i in 1:length(lines)) {
                    # check if error starts or ends
                    if(grepl("^====", lines[i])) {
                        if(!error.started) {
                            error.started <- TRUE
                        } else {
                            error.started <- FALSE
                            return(paste("Error 9: Galaxy workflow error.", error.message))
                        }
                    } else if(error.started & lines[i] != "") {
                        current.message <- paste(unlist(strsplit(lines[i], c(" ")))[-1], collapse=" ")
                        current.message <- gsub("^\\s*", " ", current.message)
                        if(error.message == "") {
                            error.message <- current.message
                        } else {
                            error.message <- paste(error.message, current.message, sep=";")
                        }
                    }
                }
            } else {
                return(paste("Error 10: Could not download error file from",url))
            }
        }
    }
    result.files <- dir(paste(factor.folder, "results", "sacurine", sep="/"))
    diff.files.indexes <- grep("\\.diff$", result.files)
    if(length(diff.files.indexes) != 4) {
        return(paste("Error 11: Unknown error. Expected 4 diff files in", paste(factor.folder, "results", "sacurine", sep="/"), ". Found only", length(diff.files.indexes)))
    }
    return("")
}

# main function does the job
main <- function (sacurine.script, study.name, log.file, study.path, wft4galaxy.template.yaml, path.to.ga.template, galaxy.key, galaxy.url, output, delete.dirs = TRUE) {
    path <- paste(study.path, study.name, sep = "/")
    # read study
    data = tryCatch({
        readISAtab(path)
    }, warning = function(w) {
        
    }, error = function(e) {
        write(paste("\"", study.name, "\",", "\"\"," , "\"\",", "\"",paste("Error 1: Could not read ISAtab folder properly.", sep=""),"\",", "\"", "", "\"", sep = ""), file=log.file,append=TRUE)
        return()
    }, finally={
        
    })
    if(is.null(data)) {
        write(paste("\"", study.name, "\",", "\"\"," , "\"\",", "\"",paste("Error 1: Could not read ISAtab folder properly.", sep=""),"\",", "\"", "", "\"", sep = ""), file=log.file,append=TRUE)
        return()
    }
    # get properties from investigation file
    study.factor.name.index <- grep("Study Factor Name", as.character(attributes(data)[["investigation.file"]][,1]))
    study.file.name.index <- grep("Study File Name", as.character(attributes(data)[["investigation.file"]][,1]))
    assay.file.name.index <- grep("Study Assay File Name", as.character(attributes(data)[["investigation.file"]][,1]))

    # get study files
    s_files<-list()
    a <- sapply(attributes(data)[["investigation.file"]][study.file.name.index,-1], function(x) {
        current.s_file <- as.character(x)
        if(!is.na(current.s_file) & current.s_file != "") {
            s_files[[length(s_files) + 1]] <<- current.s_file
        }
    })
    s_files<-unlist(s_files)
    if(length(s_files) > 1) {
        write(paste("\"", study.name, "\",", "\"", assay,"\"," , "\"",factors[factor.index],"\",", "\"",paste("Error 2: More than one s_file found", sep=""),"\",", "\"", "", "\"", sep = ""), file=log.file,append=TRUE)
        return()
    }

    # collect factors
    factors<-list()
    a <- sapply(attributes(data)[["investigation.file"]][study.factor.name.index,-1], function(x) {
        current.factor <- as.character(x)
        if(!is.na(current.factor) & current.factor != "") {
            factors[[length(factors) + 1]] <<- current.factor
        }
    })
    factors <- unlist(factors)
    if(is.null(factors)) {
        write(paste("\"", study.name, "\",", "\"\"," , "\"\",", "\"",paste("Error 2: No factors found in investigation file.", sep=""),"\",", "\"", "", "\"", sep = ""), file=log.file,append=TRUE)
        return()
    }
    # get assays
    assays<-list()
    a <- sapply(attributes(data)[["investigation.file"]][assay.file.name.index,-1], function(x) {
        current.assay <- as.character(x)
        if(!is.na(current.assay) & current.assay != "") {
            assays[[length(assays) + 1]] <<- current.assay
        }
    })
    assays <- unlist(assays)

    # get real factor names
    s_file_data <- read.csv(paste(path,s_files[1],sep="/"), sep="\t", check.names=F)
    real.factors.index <- sapply(factors, function(x) grep(paste("\\[",x,"\\]",sep=""), colnames(s_file_data)))
    
    real.factors <- list()
    a <- sapply(real.factors.index, function(x) {
        # if not found in the s file
        if(length(x) == 0) {
            real.factors[[length(real.factors) + 1]] <<- -1
        } else {
            real.factors[[length(real.factors) + 1]] <<- colnames(s_file_data)[x]
        }
    })
    real.factors <- unlist(real.factors)
    
    # get study folder for the current MTBLS study
    study.output.folder <- paste(output, study.name, sep="/")
    # remove if the current study folder already exists
    if(file.exists(study.output.folder)) {
        unlink(study.output.folder, recursive = T)
    }
    # create study folder for the current MTBLS study
    dir.create(study.output.folder)

    # run saccurine
    for(assay in assays) {
        is.valid.ms <- check.assay.file(paste(path,assay,sep="/"))
        if(!is.valid.ms) {
            write(paste("\"", study.name, "\",", "\"", assay,"\"," , "\"\",", "\"",paste("Error 4: ", assay, " is no MS study.", sep=""),"\",", "\"", "", "\"", sep = ""), file=log.file,append=TRUE)
            next
        }
        assay.folder <- paste(study.output.folder, gsub("\\.txt$","",assay), sep="/")
        dir.create(assay.folder)
        sapply(1:length(real.factors), function(factor.index) {
            success <- ""
            if(real.factors[factor.index] == "-1") {
                write(paste("\"", study.name, "\",", "\"", assay,"\"," , "\"",factors[factor.index],"\",", "\"",paste("Error 4: Did not find column in s_file: ", factors[factor.index], sep=""),"\",", "\"", "", "\"", sep = ""), file=log.file,append=TRUE)
                success <- "error"
            }
            if(success == "") {
                factor.folder <- paste(assay.folder, gsub("\\/", "_", gsub("\\s+", "_", factors[factor.index])), sep="/")
                dir.create(factor.folder)
                # define output files
                sample.file <- paste(factor.folder,"sample-output.tsv",sep="/")
                variable.file <- paste(factor.folder,"variable-output.tsv",sep="/")
                matrix.file <- paste(factor.folder,"matrix-output.tsv",sep="/")
                # run script to generate input files
                command <- 
                    paste(
                        sacurine.script," -i '",path,"'",
                        " -s '",sample.file,"'",
                        " -v '",variable.file,"'",
                        " -m '",matrix.file,"'",
                        " -f '",assay,"'",
                        " -S '",real.factors[factor.index],"'",
                        sep=""
                    )
                # run the command
                system(command)
                # check if successful
                # variable success contains empty string if 
                if(!file.exists(sample.file) || !file.exists(variable.file) || !file.exists(matrix.file)) {
                    write(paste("\"", study.name, "\",", "\"", assay,"\"," , "\"",factors[factor.index],"\",", "\"Error 5: Problem when creating Sacurine input files\",", "\"", command, "\"", sep = ""), file=log.file, append=TRUE)
                    unlink(factor.folder, recursive = T)
                    success <- "error"
                }
                factor.name <- ""
                if(success == "") {factor.name <- get.sacurine.factor.name(sample.file, factors[factor.index])}
                if(success != "" & length(factor.name) > 1) {
                    write(paste("\"", study.name, "\",", "\"", assay,"\"," , "\"",factors[factor.index],"\",", "\"",paste("Error 6: More than two matching columns found.",regex.column.name,sep=""),"\",", "\"", command, "\"", sep = ""), file=log.file, append=TRUE)
                    success <- "error"
                }
                # check for number classes
                if(success == "") {
                    success <- check.sample.file(sample.file, factor.name)
                    if(success != "") {
                        write(paste("\"", study.name, "\",", "\"", assay,"\"," , "\"",factors[factor.index],"\",", "\"",success,"\",", "\"", command, "\"", sep = ""), file=log.file, append=TRUE)
                        unlink(factor.folder, recursive = T)
                    }
                }
                # prepare, run and validate wft4galaxy for sacurine workflow
                if(success == "") {
                    prepare.wft4galaxy.files(factor.name, factor.folder, wft4galaxy.template.yaml, path.to.ga.template)
                    run.wft4galaxy(factor.folder, galaxy.key, galaxy.url)
                    success <- validate.wft4galaxy.run(factor.folder, galaxy.url)
                    if(success != "") {
                        write(paste("\"", study.name, "\",", "\"", assay,"\"," , "\"",factors[factor.index],"\",", "\"",success,"\",", "\"", command, "\"", sep = ""), file=log.file, append=TRUE)
                        # better to not delete the current folder as the wft4galaxy files might be interesting for further analysis
                    #    unlink(factor.folder, recursive = T)
                    }
                }
                # write positive feedback to log file
                if(success == "") {
                    write(paste("\"", study.name, "\",", "\"", assay,"\"," , "\"",factors[factor.index],"\",", "\"OK!\",", "\"", command,"\"", sep = ""), file=log.file, append=TRUE)
                }
            }
        })
        # remove assay folder if empty
        if(length(dir(assay.folder)) == 0) {
            if(delete.dirs) {unlink(assay.folder, recursive = T)}
        }
    }
    # remove study folder if empty
    if(length(dir(study.output.folder)) == 0) {
        if(delete.dirs) {unlink(study.output.folder, recursive = T)}
    }
}

# run the command

main(sacurine.script, study.name, log.file, study.path, wft4galaxy.template.yaml, path.to.ga.template, galaxy.key, galaxy.url, output, FALSE)
