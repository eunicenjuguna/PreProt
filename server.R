# install.packages(c("shiny","tidyverse","DT","shinydashboard","shinydashboardPlus"))
# 
# install.packages(c("tidyverse","data.table","janitor","data.table",
#                    "ggfortify", "scales","ggplot2","reshape2",
#                    "ggpubr","visdat","plotly", "dlookr","missRanger",
#                    "pcaMethods","imputeLCMD","impute","sva",
#                    "variancePartition", "pvca")) 

#   if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#   BiocManager::install(version = "3.15")
# 
#   # BiocManager::install(c("impute","sva","limma","variancePartition"))
# BiocManager::install("vsn")

library(shiny)
library(tidyverse)
library(data.table)
library(DT)
library(janitor)# load the package to easily sum columns and rows
library(data.table)
library(ggfortify)
library(scales)
library(ggplot2)
library(reshape2) 
library(ggpubr)
library(visdat)
library(plotly)
library(dlookr)
library(missRanger)
library(missForest)
library(pcaMethods)
library(imputeLCMD)
library(impute)
library(sva)
library(limma)
library(vsn)
library(variancePartition)


#vis <- read_csv("vis.csv")

server <- function(input, output) {
  options(shiny.maxRequestSize=30*1024^200)
  
  tmt_file <- reactive({
    
    #condition to input data
    req(input$tmt_file)
    ## assigning input to a object
    infile = input$tmt_file
    # actual reading
    df = read.table(infile$datapath, header = T, sep = "\t")
    df
  })
  
  output$tmt_raw_data <- renderDataTable({
    
    tmt_file() 
  }, options = list(
    
    scrollX = TRUE))
  
  #autoWidth = F,)
  
  tmt_corrected_intensity_file <- reactive({
    
    reporter_intesinty <-  tmt_file() %>%
      mutate(protein.id = gsub(";.*$","", Protein.IDs)) %>%
     mutate(protein.name = gsub(";.*$","", Protein.names)) %>%  #select from;and replace with empty string
      select(protein.id,protein.name,
             #matches("Reporter\\.intensity\\.corrected\\.[0-9]{1,2}\\.TMT[0-9]{1,2}$"))
             matches("Reporter\\.intensity\\.corrected\\.[0-9]{1,2}\\.[0-9]{1,2}$|Reporter\\.intensity\\.corrected\\.[0-9]{1,2}\\.[[:alpha:]]{1,9}[0-9]{1,2}$"))#use it if you require reporter intensity corrected
             #matches("Reporter\\.intensity\\.[0-9]"))
   #  tmt_col= names(reporter_intesinty)
   # old_tmt= tmt_col[str_detect(tmt_col,"Reporter\\.intensity\\.corrected")]
   # item1=old_tmt[1]
   # if(str_detect(item1,"Reporter\\.intensity\\.corrected\\.[0-9]{1,2}\\.[[:alpha:]]{1,9}[0-9]{1,2}$")){
   #   batch_no= str_extract(old_tmt,"[0-9]{1,2}$")
   #   new_tmt= gsub("\\.[0-9]{1,2}\\.[[:alpha:]]{1,9}[0-9]{1,2}$","",old_tmt)
   #   new_tmt= paste0(new_tmt,".",batch_no)
   #   setnames(reporter_intesinty,old_tmt,new_tmt)
   #}

    reporter_intesinty
    
  })
  

  tmt_reporter_intensity_file <- reactive({
    
    reporter_intesinty <-  tmt_file() %>%
      mutate(protein.id = gsub(";.*$","", Protein.IDs)) %>% 
      mutate(protein.name = gsub(";.*$","", Protein.names)) %>% 
      select(protein.id,protein.name,
             matches("Reporter\\.intensity\\.[0-9]{1,2}\\.[0-9]{1,2}$|Reporter\\.intensity\\.[0-9]{1,2}\\.[[:alpha:]]{1,9}[0-9]{1,2}$"))
            # matches("Reporter\\.intensity\\.[0-9]{1,2}\\.TMT[0-9]{1,2}$"))
             #matches("Reporter\\.intensity\\.[0-9]{1,2}\\.[0-9]{1,2}$"))
    
    reporter_intesinty
    
  })
  
  ##Choose corrected or uncorrected intensities
  # 
  choose_intesity <- reactive({
    if(input$reporter_intesi_corrected == "Reporter Intensity Corrected"){

      tmt_corrected_intensity_file()
    }else{
      tmt_reporter_intensity_file()
    }
  })

  output$corrected_title <- renderText({

    if(input$reporter_intesi_corrected == "Reporter Intensity Corrected"){

      intesity_title = "Reporter intesinty corrected"
    }else{
      intesity_title = "Reporter intesinty"
    }

    intesity_title

  })
  
  output$tmt_corrected_intensity <- renderDataTable({
    
    choose_intesity() 
  }, options = list(
    
    scrollX = TRUE))
  
  output$down_corrected_intensity_tmt <- downloadHandler(
    filename = function() {
      paste0("reporter_intenstity_info",".csv")
    },
    content = function(file) {
      fwrite(choose_intesity(), file, row.names = FALSE)
    }
  )
  
  

  
  
  lfq_file <- reactive({
    
    #condition to input data
    req(input$lfq_file)
    ## assigning input to a object
    infile = input$lfq_file
    # actual reading
    df = read.table(infile$datapath, header = T, sep = "\t")
    df
  })
  
  output$lfq_raw_data <- renderDataTable({
    
    lfq_file() 
  }, options = list(
    
    scrollX = TRUE))
  
  #autoWidth = F,)
  
  lfq_intensity_data_file <- reactive({
    
    reporter_intesinty <-  lfq_file()%>%
      mutate(protein.id = gsub(";.*$","", Protein.IDs)) %>% #select from;and replace with empty
      mutate(protein.name = gsub(";.*$","", Protein.names)) %>% 
      select(protein.id,protein.name,
             contains("LFQ.intensity"))
    
    reporter_intesinty
    
  })
  
  
  output$lfq_intensity_data <- renderDataTable({
    
    lfq_intensity_data_file() 
  }, options = list(
    
    scrollX = TRUE))
  
  output$down_lfq_intensity_lfq <- downloadHandler(
    filename = function() {
      paste0("lfq_intensity_info",".csv")
    },
    content = function(file) {
      fwrite(lfq_intensity_data_file(), file, row.names = FALSE)
    }
  )
  
  
choose_analysis_type<- reactive({
  if(input$select_analysis_type == "TMT" ){
    choose_intesity()
    
  }else{
    lfq_intensity_data_file() 
  }
  
}) 
  
  
  threshold_data_file <- reactive({
    reporter_intesinty <- choose_analysis_type()
    number_of_zeros <- rowSums( reporter_intesinty==0)
    percentage_missing <- round(number_of_zeros/ncol(reporter_intesinty) * 100, 2)
    
    reporter_intesinty <- reporter_intesinty  %>%
      #ProteinGroups <- ProteinGroups %>%
      mutate(percentage_missing = percentage_missing)
    
    ## remove rows (observations)with more than 80% missing data
    ProteinGroups_data_sub <- reporter_intesinty %>%
      filter(percentage_missing < input$filter_threshold)
    #row1 <- reporter_intesinty[1,] %>%  as.numeric()
    
    ProteinGroups_data_sub
  })
    
    output$threshold_data <- renderDataTable({
      
      threshold_data_file() 
    }, options = list(
      
      scrollX = TRUE))
    
    output$down_threshold_data <- downloadHandler(
      filename = function() {
        paste0("threshold_filtered",".csv")
      },
      content = function(file) {
        fwrite(threshold_data_file() , file, row.names = FALSE)
      }
    )
    
    
    cont_free_data_file <- reactive({
    #Remove all the observations containing the contaminants (keratin,amalyse,reverse)
      threshold_data =  threshold_data_file()
      threshold_data<- threshold_data %>%
        select(-percentage_missing)
      contaminats =input$filter_contaminants %>% tolower()
      contaminats_p <- paste0(contaminats, collapse = "|")
      without_last <- threshold_data %>%
        filter(!grepl(contaminats_p, tolower(protein.name))) %>%
        filter(!grepl(contaminats_p, tolower(protein.id)))
        
    # 
    # #cut all observations containing REV (reverse sequences) and !dont show them
    # without_REV <-  without_Keratin %>%
    #   filter(!grepl("REV", protein.id))
    # 
    # #cut all observations containing CON (contaminats)and !dont show them
    # without_CON <-  without_REV %>%
    #   filter(!grepl("CON", protein.id))
    # 
    # without_amalyse <-  without_CON %>%
    #   filter(!grepl("amalyse", protein.id))
    # 
    # 
    # without_last<- without_CON %>%
    #   select(-percentage_missing)
    
    without_last
   
    
  })
  
  output$cont_free_data <- renderDataTable({
    
    cont_free_data_file() 
  }, options = list(
    
    scrollX = TRUE))
  
  output$down_cont_free_data <- downloadHandler(
    filename = function() {
      paste0("contaminants_free",".csv")
    },
    content = function(file) {
      fwrite( cont_free_data_file() , file, row.names = FALSE)
    }
  )
  

  missing_values_inspection_file <- reactive({
    meta_data=cont_free_data_file()[,c(1:3)]
    intensity_data = cont_free_data_file()[,-c(1:3)]
    intensity_data[intensity_data == 0] <- NA
    #inspect the missingness
    
    intensity_data
   
  })
  mv_inspect<- reactive({
    inspect_mv<- plot_na_pareto(missing_values_inspection_file(),only_na = TRUE,
                                grade = list(low = 0.2, high = 0.5,  very_high = 1),
                                main = "The Percentage of Missing Values in your Data",plot = FALSE)
    inspect_mv %>%
      mutate(percentage= round(ratio*100, 2))%>%
      select(-cumulative,-ratio)%>%
      arrange(desc(percentage))
    
  })
  output$missing_values_inspection <- renderDataTable({
    mv_inspect()
    
    
   
   })
  output$down_inspect_mv <- downloadHandler(
    filename = function() {
      paste("missing_values_percentage", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(mv_inspect(), file, row.names = FALSE)
    })
 
  
  missing_values_type_file <- reactive({
    
    model_select <-model.Selector( missing_values_inspection_file())
    model_select 
    
  })
  
  output$missing_values_type<-renderText({
   model_select<- missing_values_type_file()
   
   mv_type<- model_select[[1]]%>%
     sum()
  mv_row<- missing_values_inspection_file()%>%
     nrow()
  mnar_calc<- mv_type/ mv_row * 100
  if(mnar_calc >= 80 ){
   
    "<p> <strong> The missing values in your data are missing at random(MAR)and missing
    completely at random(MCAR). Select your preferred imputation method in the 
    next step from the first three methods namely; KNN, missRanger or missForest.</strong>
    <br>
        This step serves to inform the user about the type of missing values in the
    data. The data is expected to be missing at random (MAR), missing completely
    at random (MCAR), or missing not at random (MNAR). Proteomics studies suggest
    that the type of missing values in data determines the imputation method to
  be applied (Donders et al., 2006; Gardner & Freitas, 2021; Wei et al., 2018).
 <br>
 PreProt employed this approach by informing the user on the imputation method
  to use depending on the type of missing values in a given dataset.The next step
  is to impute the missing values using four imputation methods incorporated 
  in PreProt. K nearest neighbor (KNN) which is the default method,
two random forest (RF) methods namely missRanger and missForest and finally 
quantile regression imputation of left censored data (QRILC) method. KNN, 
missRanger and missForest imputation methods are known to work well with data
with missing values missing at random and completely at random while QRILC is
known to impute data missing not at random. Therefore, the user is advised to 
select an imputation method depending on the type of missingness in the data.
    </p>"


    
  }else{
    "<p> <strong> The missing values in your data are missing not at random(MNAR).
    Select QRLIC imputation method imputation method in the next step.</strong> 
        <br>
          This step serves to inform the user about the type of missing values in the
        data. The data is expected to be missing at random (MAR), missing completely
        at random (MCAR) or missing not at random (MNAR). Proteomics studies suggest
        that the type of missing values in data determines the imputation method to
        be applied (Donders et al., 2006; Gardner & Freitas, 2021; Wei et al., 2018).
        <br>
          PreProt employed this approach by informing the user on the imputation method
        to use depending on the type of missing values in a given data.The next step
        is to impute the missing values using four imputation methods incorporated 
        in PreProt. K nearest neighbor (KNN) which is the default method,
        two random forest (RF) methods namely missRanger and missForest and finally 
        quantile regression imputation of left censored data (QRILC) method. KNN, 
        missRanger and missForest imputation methods are known to work well with data
        with missing values missing at random and completely at random while QRILC is
        known to impute data missing not at random. Therefore, the user is advised to 
        select an imputation method depending on the type of missingness in the data.
        </p>"

    
  }
  }
   
  )

  missRanger_file <- reactive({
    imputed_ranger<-
      missRanger(
        missing_values_inspection_file(),
        formula = . ~ .,
        pmm.k = 0L,
        maxiter = 10L,
        seed = 3,
        num.tree = 1000,
      )
    meta_data=cont_free_data_file()[,c(1:3)]
    imputed_data_ranger = cbind(meta_data[,c(1:3)], imputed_ranger)
    
    imputed_data_ranger
    
  })
  
  missForest_file <- reactive({
    meta_data=cont_free_data_file()[,c(1:3)]
    imputed_forest<- missForest(xmis = missing_values_inspection_file(), maxite= 100,ntree = 1000, maxnodes = 2)
    
    imputed_forest<-imputed_forest$ximp
    
    imputed_data_forest = cbind(meta_data[,c(1:3)],imputed_forest )
    
    imputed_data_forest
    
  })
  KNN_file <- reactive({
    meta_data=cont_free_data_file()[,c(1:3)]
    imputeby_knn <- impute.knn(as.matrix(missing_values_inspection_file()), k = 10, rowmax = 0.5,colmax = 1, 
                               maxp = 1500, rng.seed=362436069)
    
    imputed_knn<- imputeby_knn$data 
    
    #intensity_data = ProteinGroups_clean[,-c(1:4)]
    #meta_data= ProteinGroups_clean[,c(1:4)]
    imputed_data_knn = cbind(meta_data[,c(1:3)], imputed_knn)%>%
      as.data.frame()
    
    imputed_data_knn 
    
  })
  
  qrilc_file <- reactive({
    meta_data=cont_free_data_file()[,c(1:3)]
    qrilc_data<- impute.QRILC(missing_values_inspection_file(), tune.sigma = 1)
    qrilc<- qrilc_data[[1]]
  
    
    #intensity_data = ProteinGroups_clean[,-c(1:4)]
    #meta_data= ProteinGroups_clean[,c(1:4)]
    imputed_data_qrilc = cbind(meta_data[,c(1:3)],qrilc )%>%
      as.data.frame()
    
    imputed_data_qrilc 
    
  })
  
  
  
  choose_imputation_method_file <- reactive({
    if(input$choose_imputation_method == "missRanger"){
      
      missRanger_file()
    }else if(input$choose_imputation_method == "KNN"){
      KNN_file ()
    }
    else if(input$choose_imputation_method == "QRILC"){
      qrilc_file ()
    }else{
      missForest_file()
    }
  })
  
  output$imputation_title <- renderText({
    
    if(input$choose_imputation_method == "missRanger_file"){
      
      imputation_title = "Machine learning missRanger imputation method"
    }else if(input$choose_imputation_method == "KNN"){
      imputation_title = " K nearest neigbour (KNN) imputation method"
    }
    else if(input$choose_imputation_method == "QRILC"){
      imputation_title = " QRILC imputation method for MNAR Values"
    }else{
       imputation_title = "Machine learning missForest imputation method "
    }
    
    imputation_title
    
  })
  
  output$imputation_method<- renderDataTable({
    
    
    choose_imputation_method_file()
  }, options = list(
    
    scrollX = TRUE))
  
  output$down_down_imputation_method <- downloadHandler(
    filename = function() {
      paste0("imputation_results",".csv")
    },
    content = function(file) {
      fwrite( choose_imputation_method_file (), file, row.names = FALSE)
    }
  )
  
  
  
  before_norm_file <- reactive({
    meta_data=cont_free_data_file()[,c(1:2)]
   # intensity_data = choose_imputation_method_file ()[,-c(1:3)]
    #intensity_data[intensity_data == -Inf] <- NA
   before_norm_data<- choose_imputation_method_file ()[,-c(1,2)]
   before_norm_data
  })
  output$before_norm <- renderPlot({
    boxplot(before_norm_file())
    
  })
 
    
  log_normalization_file <- reactive({
    
    meta_data=cont_free_data_file()[,c(1:2)]
    #intensity_data = cont_free_data_file()[,-c(1:3)]
    #intensity_data[intensity_data == -Inf] <- NA
    
    log_norm<- log(choose_imputation_method_file ()[,-c(1,2)])
    
    log_normalized<- cbind(meta_data[,c(1:2)], log_norm)
    
    log_normalized
    
  })
  
  
  vsn_normalization_file <- reactive({
    
    meta_data=cont_free_data_file()[,c(1:2)]
    # intensity_data = cont_free_data_file()[,-c(1:3)]
    # intensity_data[intensity_data == -Inf] <- NA
    # 
    vsn_norm = normalizeVSN(choose_imputation_method_file ()[,-c(1,2)])
    
    vsn_normalized<- cbind(meta_data[,c(1:2)], vsn_norm)
    
    vsn_normalized
  })

  cyclicloess_file<- reactive({
    
    meta_data=cont_free_data_file()[,c(1:2)]
    # intensity_data = cont_free_data_file()[,-c(1:3)]
    # intensity_data[intensity_data == -Inf] <- NA
    # 
    log_norm<- log(choose_imputation_method_file ()[,-c(1,2)])
    
    loess_norm = normalizeCyclicLoess(log_norm)
    
    loess_normalized<- cbind(meta_data[,c(1:2)], loess_norm)
    
    loess_normalized
  })
  quantile_file<- reactive({
    qant<- choose_imputation_method_file ()[,-c(1,2)]
    meta_data=cont_free_data_file()[,c(1:2)]
    
    quantile_norm = normalizeQuantiles(qant)
    quantile_normalized<- cbind(meta_data[,c(1:2)], quantile_norm)
    
    quantile_normalized
  })
  
  mean_file<- reactive({
    
    meta_data=cont_free_data_file()[,c(1:2)]
    # intensity_data = cont_free_data_file()[,-c(1:3)]
    # intensity_data[intensity_data == -Inf] <- NA
    mean_norm_minmax <- function(x){
      (x- mean(x)) /(max(x)-min(x))
    }
    mean_norm <- as.data.frame(lapply(choose_imputation_method_file ()[,-c(1,2)], mean_norm_minmax))
    #mean_norm = normalizeQuantiles(choose_imputation_method_file ()[,-c(1,2)])
    mean_normalized<- cbind(meta_data[,c(1:2)], mean_norm)
    
    mean_normalized
  }) 
  ##Choose type
  # 
  choose_normalization <- reactive({
    if(input$normalization_method == "Log Normalization"){
      
      
      log_normalization_file()
    }else if(input$normalization_method == "VSN"){
      vsn_normalization_file()
    }else if(input$normalization_method == "Cyclicloess"){
      cyclicloess_file()
    }else if(input$normalization_method == "Quantile"){
      quantile_file()
    }else{mean_file()}
  })
  
  output$normalization_title <- renderText({
    
    if(input$normalization_method == "Log Normalization"){
      
      log_title = "Log Normalization"
    }else if(input$normalization_method == "VSN"){
      log_title = "VSN Normalization"
    }else if(input$normalization_method == "Cyclicloess"){
      log_title = " Cyclicloess Normalization"
    }else if(input$normalization_method == "Quantile"){
      log_title = " Quantile Normalization"
    }else{
      log_title = " Mean Normalization"}
    
    log_title
  })
  
  output$norm_method <- renderDataTable({
    
    choose_normalization() 
  }, options = list(
    
    scrollX = TRUE))
  
  output$down_norm_method <- downloadHandler(
    filename = function() {
      paste0("log_normalizatio_results",".csv")
    },
    content = function(file) {
      fwrite(choose_normalization(), file, row.names = FALSE)
    }
  )
  after_norm_file <- reactive({
   meta_data=cont_free_data_file()[,c(1:2)]
    # intensity_data = choose_imputation_method_file ()[,-c(1:3)]
    #intensity_data[intensity_data == -Inf] <- NA
    
 
    after_norm_data<- choose_normalization()[,-c(1,2)]
    after_norm_data
  })
  output$after_norm <- renderPlot({
    boxplot(after_norm_file())
    
    
  })
  
  

  
  metadata_ouput <- reactive({
    
    #condition to input data
    req(input$metadata_file)
    ## assigning input to a object
    infile = input$metadata_file
    # actual reading
    df = read.csv(infile$datapath, header = T)
    df
    
  })
  
  output$metadata_data_df <- renderDataTable({
    
    metadata_ouput() 
  }, options = list(
    
    scrollX = TRUE))
  
  
  variation_source_file <- reactive({
    imp_exp<- choose_imputation_method_file()[-c(2)]%>%
      column_to_rownames("protein.id")%>%
      select(matches("^reporter.intensity|^Reporter.intensity"))
    
    imp_exp_t <- t(imp_exp)
    
    # info <- read.csv("tmt information.csv", header = T)
    match_rownames<- rownames(imp_exp_t)
    info_info <- cbind(match_rownames,metadata_ouput() )
    pvca_batch_info <- info_info[,-1]
    rownames(pvca_batch_info) <- info_info[,1]
    
    
    
    
    # pvca_info <- cbind(imp_exp_t,info)
    
    all(rownames(pvca_batch_info) == colnames(imp_exp))
    
    
    form <- ~(1|batch) + (1|run_day)+ (1|sex) + (1|case_control) + (1|sample_source) +(1|sample_type) +  (1|sampling_date) 
    
    pvca_batch_info$batch<- as.factor(pvca_batch_info$batch)
    pvca_batch_info$run_day<- as.factor(pvca_batch_info$run_day)
    pvca_batch_info$sex<- as.factor(pvca_batch_info$sex)
    pvca_batch_info$case_control<- as.factor(pvca_batch_info$case_control)
    pvca_batch_info$sample_source<- as.factor(pvca_batch_info$sample_source)
    pvca_batch_info$sample_type<- as.factor(pvca_batch_info$sample_type)
    pvca_batch_info$sampling_date<- as.factor(pvca_batch_info$sampling_date)
    
    
    
    
    # browser()
    varPart <- fitExtractVarPartModel( imp_exp, form, pvca_batch_info)
    
    vp <- sortCols( varPart )
    plotPercentBars( vp[1:10,] )
    # plotVarPart( vp )
    
    vp
    
  })
  output$variation_source_plot <- renderPlot({
    plotVarPart(variation_source_file())
    
    
  })
  
  
  inspect_batcheffect_file <- reactive({
    
    library(data.table)
    # validate(
    #need(try(input$select_analysis_type == "TMT"), "Batch Correction Applies for TMT Data")
    # )
    
    # req(input$select_analysis_type == "TMT")
    
    inspect_batch <- choose_normalization()
    
    reporter_tras <- transpose(inspect_batch,keep.names = "intensity")#tranpose it for plotting
    reporter_id_nms <- reporter_tras  %>% 
      filter(intensity  == "protein.id") %>% #select what will be the header
      select(-intensity)#do not filter the intesity
    reporter_names <- names(reporter_id_nms)#make the ids the heads
    reporter_new_names <- reporter_id_nms[1,] %>%#make them characters
      as.character()
    
    
    #set the new names (protein ids) to be headers
    setnames(reporter_tras , reporter_names, reporter_new_names)#set the names from old to new
    
    reporter_tras  <- reporter_tras  %>% #in the tranposed observation 
      # create onother column ((batch_sample)remember reporter intensity corrected columns 
      #contained numbers i.e(sample and batcg) which will bring trouble plotting the graph
      mutate(sample_batch = str_extract(intensity, "[0-9]{1,2}_[0-9]{1,2}|[0-9]{1,2}\\.[0-9]{1,2}$|[0-9]{1,2}\\.[[:alpha:]]{1,9}[0-9]{1,2}$"))%>%
      mutate(sample_batch=gsub("[[:alpha:]]{1,9}","",sample_batch))
    
    reporter_tras  <- reporter_tras  %>% #
      filter(!intensity %in%
               c("protein.name","protein.id"))
    reporter_tras  <- reporter_tras %>% 
      mutate(sample_batch= gsub("_","\\.",sample_batch)) %>%
      separate(sample_batch, into = c("sample", "batch"), sep = "\\.",
               remove = F)# now separate
    #remove the row x
    reporter_tras  <- reporter_tras %>%
      filter(intensity!= "X.1")
    
    reporter_tras  <- reporter_tras %>%
      filter(intensity!= "X")
    
    reporter_tras <-  reporter_tras %>%
      select(intensity,sample_batch,sample,batch,all_of(reporter_new_names))
    
    reporter_tras_batch<- reporter_tras %>%
      select(sample_batch,sample,batch)
    
    
    #tmt_info <- read.csv("tmt information.csv", header = T)
    
    # select only the colums of interest,(the proteins) 
    reporter_tras_pca  <- reporter_tras %>%
      select(all_of(reporter_new_names))
    
    
    reporter_tras_pca<- reporter_tras_pca %>%
      mutate_at(reporter_new_names, as.numeric)
    
    pca_data <- prcomp(reporter_tras_pca, scale. = TRUE)
    
    #autoplot(pca_data)
    pca_10 <- cbind(reporter_tras_batch,pca_data$x[,1:5])#create pc 1-10
    #plot PC1 and PC2
    #library("ggplot2") # load the required plotting package
    
    #pca_10 <- pca_10%>%
    # filter(sample!=10)
    pca_sdev <- pca_data$sdev^2
    
    pca_sdev_per <- round(pca_sdev/sum(pca_sdev)*100, 1)
    
    pca_10
    
    list(pca_10 = pca_10,
         pca_sdev_per = pca_sdev_per,
         pca_sdev = pca_sdev,
         reporter_tras = reporter_tras,
         reporter_tras_batch = reporter_tras_batch)
    
  })
  output$inspect_batcheffect <- renderPlot({
    pca_10 = inspect_batcheffect_file()$pca_10%>%
      mutate(batch=as.factor(batch))
    pca_sdev_per = inspect_batcheffect_file()$pca_sdev_per
    ggplot(data = inspect_batcheffect_file()$pca_10, 
           aes(x = PC1, y = PC2 , color = batch, 
               group = batch)) +
      #aes(x = PC1, y = PC2 , color = as.factor(run_day), group = as.factor(run_day))) +
      geom_point() +   geom_point(size=1)+
      theme_pubr(border = TRUE,legend = 'right')+
      theme(text = element_text(size = 30))+
      
      xlab(paste("PC1 - ", pca_sdev_per[1], "%", sep="")) +
      
      ylab(paste("PC2 - ", pca_sdev_per[2], "%", sep="")) +
      
      theme_bw() +
      
      ggtitle("PCA before batch correction")+
      
      stat_ellipse(type = "norm")
    
  })
  
  
  batch_correction_file <- reactive({
    
    meta_data=cont_free_data_file()[,c(1:2)]
    reporter_tras_batch = inspect_batcheffect_file()$reporter_tras_batch
    
    batch = as.factor(reporter_tras_batch$batch)
    # # #modcombat = model.matrix(~1,data = imputed_ranger)
    combat_n = ComBat(dat = choose_normalization()[,-c(1,2)], batch=batch,
                      par.prior = TRUE, prior.plots = FALSE)
    
    corrected<- cbind(meta_data[,c(1:2)], combat_n)
    #corrected<-  choose_normalization()[,-c(1,2)]
    
    corrected
    
  })
  output$batch_correction <- renderDataTable({
    
    batch_correction_file () 
  }, options = list(
    
    scrollX = TRUE))
  output$down_batch_correction <- downloadHandler(
    filename = function() {
      paste0("batch corrected",".csv")
    },
    content = function(file) {
      fwrite(  batch_correction_file (), file, row.names = FALSE)
    }
  )
  
  visualize_after_correction_file <- reactive({
    reporter_tras<- inspect_batcheffect_file()$reporter_tras
    sample_batch<- reporter_tras[,1:4]
    corrected_tras <- transpose(batch_correction_file(),keep.names = "intensity")#tranpose it for plotting
    
    corrected_id_nms <- corrected_tras  %>% 
      filter(intensity  == "protein.id") %>% 
      select(-intensity)#do not filter the intesity
    corrected_names <- names(corrected_id_nms)#make the ids the heads
    corrected_new_names <- corrected_id_nms[1,] %>%#coarse characters
      as.character()
    
    setnames(corrected_tras , corrected_names, corrected_new_names)#set the names from old to new
    
    corrected_tras  <- corrected_tras  %>% 
      mutate(sample_batch = str_extract(intensity,"[0-9]{1,2}_[0-9]{1,2}|[0-9]{1,2}\\.[0-9]{1,2}$|[0-9]{1,2}\\.[[:alpha:]]{1,9}[0-9]{1,2}$"))%>%
      mutate(sample_batch=gsub("[[:alpha:]]{1,9}","",sample_batch))
    # mutate(sample_batch = str_extract(intensity, "[0-9]{1,2}\\.[0-9]{1,2}$|[0-9]{1,2}\\.[[:alpha:]]{1,9}[0-9]{1,2}$"))%>%
    # mutate(sample_batch=gsub("[[:alpha:]]{1,9}","",sample_batch))
    
    #mutate(sample_batch = str_extract(intensity, "[0-9]{1,2}\\.[0-9]{1,2}$"))
    
    corrected_tras  <- corrected_tras  %>% 
      filter(!intensity %in%
               c("Protein.name","protein.id","protein.name"))
    corrected_tras  <- corrected_tras %>% 
      mutate(sample_batch = gsub("_", ".", sample_batch))%>%
      separate(sample_batch, into = c("sample", "batch"), sep = "\\.",
               remove = F)
    corrected_tras  <- corrected_tras %>%
      filter(intensity!= "X.1")
    
    corrected_tras  <- corrected_tras %>%
      filter(intensity!= "X")
    
    corrected_tras_batch<- corrected_tras %>%
      select(sample_batch,sample,batch)
    
    # select only the colums of interest,(the proteins) 
    corrected_tras_pca  <- corrected_tras %>%
      select(all_of(corrected_new_names))
    
    # write.csv(corrected_tras, "corr_trial.csv")    
    corrected_tras_pca<- corrected_tras_pca %>%
      mutate_at(corrected_new_names, as.numeric)
    pca_data2 <- prcomp(corrected_tras_pca)
    
    #autoplot(pca_data)
    pca_10_c <- cbind(corrected_tras_batch,pca_data2$x[,1:5])#create pc 1-10
    
    #plot PC1 and PC2
    #library("ggplot2") # load the required plotting package
    
    #pca_10 <- pca_10%>%
    # filter(sample!=10)
    
    pca_sdev <- pca_data2$sdev^2
    
    pca_sdev_per <- round(pca_sdev/sum(pca_sdev)*100, 1)
    
    
    list(pca_sdev=pca_sdev,
         pca_sdev_per=pca_sdev_per,
         pca_10_c=pca_10_c)
  })
  output$visualize_after_correction <- renderPlot({
    
    pca_10_c=visualize_after_correction_file ()$pca_10_c%>%
      mutate(batch=as.factor(batch))
    # write_csv(pca_10_c, file = "LFQ afterB.csv")
    pca_sdev_per=visualize_after_correction_file ()$pca_sdev_per
    ggplot(data = pca_10_c,
           #aes(x = PC1, y = PC2 , color = as.factor(sample), group = sex)) +
           aes(x = PC1, y = PC2 , color = batch
               , group = batch)) +
      geom_point() + 
      
      xlab(paste("PC1 - ", pca_sdev_per[1], "%", sep="")) +
      
      ylab(paste("PC2 - ", pca_sdev_per[2], "%", sep="")) +
      
      theme_bw() +
      
      ggtitle("PCA after batch correction")+
      
      stat_ellipse(type = "norm")
  })
  
  other_variation_sources_file <- reactive({
    corrected_exp <-  batch_correction_file()[-c(2)]%>%
      column_to_rownames("protein.id")
    # 
    # corrected_exp_final<- corrected_exp%>%
    #   column_to_rownames("protein.id")%>%
    #   select(matches("^reporter.intensity|^Reporter.intensity"))
    
    imp_exp_t <- t(corrected_exp)
    
    
    match_rownames<- rownames(imp_exp_t)
    info_info <- cbind(match_rownames,metadata_ouput() )
    pvca_batch_info <- info_info[,-1]
    rownames(pvca_batch_info) <- info_info[,1]
    
    
    
    all(rownames(pvca_batch_info) == colnames(corrected_exp))
    
    
    form <- ~(1|batch) + (1|run_day)+ (1|sex) + (1|case_control) + (1|sample_source) +(1|sample_type) +  (1|sampling_date) 
    
    pvca_batch_info$batch<- as.factor(pvca_batch_info$batch)
    pvca_batch_info$run_day<- as.factor(pvca_batch_info$run_day)
    pvca_batch_info$sex<- as.factor(pvca_batch_info$sex)
    pvca_batch_info$case_control<- as.factor(pvca_batch_info$case_control)
    pvca_batch_info$sample_source<- as.factor(pvca_batch_info$sample_source)
    pvca_batch_info$sample_type<- as.factor(pvca_batch_info$sample_type)
    pvca_batch_info$sampling_date<- as.factor(pvca_batch_info$sampling_date)
    
    
    varPart <- fitExtractVarPartModel( corrected_exp, form, pvca_batch_info)
    
    
    vp <- sortCols( varPart )
    plotPercentBars( vp[1:10,] )
    
    vp
    
  })
  output$other_variation_source <- renderPlot({
    plotVarPart( other_variation_sources_file())
    
    
  })
  
  
}
