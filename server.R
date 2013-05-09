# shiny application # runApp("~/Documents/shinyApp")
library(shiny)
# violin plot
library(vioplot)
# density plot
library(lattice)


# Mutation visualization
mutChar = c(21,23,24,25,22)
mutCol = c(500,30,500,30,500)
mutName = matrix(NA,ncol=2,nrow=7)
mutName[,1] = c("Missense_Mutation","Nonsense_Mutation","Splice_Site","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins")
mutName[,2] = c(1,2,3,4,4,5,5)

# RAS - Expression and copy number
ras.exp = read.csv("/Users/senbabay/Documents/shinyApp/data/RSEM_Chr22_RAS.txt",header=T,sep="\t")
ras.cn = read.csv("/Users/senbabay/Documents/shinyApp/data/CN_GeneBySample_Chr22_RAS.txt",header=T,sep="\t")
ras.mut = read.csv("/Users/senbabay/Documents/shinyApp/data/Mut_Chr22_RAS.txt",header=T,sep="\t")

mx = match(substr(colnames(ras.exp),1,12),substr(colnames(ras.cn),1,12))
ras.cn2 = ras.cn[,mx[!is.na(mx)]]
ras.exp2 = ras.exp[,which(!is.na(mx))]

# BRAF - Expression and copy number
braf.exp = read.csv("/Users/senbabay/Documents/shinyApp/data/RSEM_Chr22_BRAF.txt",header=T,sep="\t")
braf.cn = read.csv("/Users/senbabay/Documents/shinyApp/data/CN_GeneBySample_Chr22_BRAF.txt",header=T,sep="\t")
braf.mut = read.csv("/Users/senbabay/Documents/shinyApp/data/Mut_Chr22_BRAF.txt",header=T,sep="\t")

mx = match(substr(colnames(braf.exp),1,12),substr(colnames(braf.cn),1,12))
braf.cn2 = braf.cn[,mx[!is.na(mx)]]
braf.exp2 = braf.exp[,which(!is.na(mx))]
rm(braf.cn)
rm(braf.exp)

# Methylation file - BRAF
braf.meth = read.csv("/Users/senbabay/Documents/shinyApp/data/chr22_bestcor_probes_sample_data_BRAF.txt",header=T,sep="\t")
mxb = match(substr(colnames(braf.exp2),1,12),colnames(braf.meth)[7:ncol(braf.meth)]) +6
braf.meth2 = braf.meth[,mxb[!is.na(mxb)]]
rownames(braf.meth2) = braf.meth[,4]
rm(braf.meth)

# Methylation file - BRAF
ras.meth = read.csv("/Users/senbabay/Documents/shinyApp/data/chr22_bestcor_probes_sample_data_RAS.txt",header=T,sep="\t")
mxr = match(substr(colnames(ras.exp2),1,12),colnames(ras.meth)[7:ncol(ras.meth)]) +6
ras.meth2 = ras.meth[,mxr[!is.na(mxr)]]
rownames(ras.meth2) = ras.meth[,4]
rm(ras.meth)

# Master file
masterB  = read.csv("/Users/senbabay/Documents/shinyApp/data/masterResultsMatrix_BRAF_v2.txt",sep="\t",header=TRUE)
masterR  = read.csv("/Users/senbabay/Documents/shinyApp/data/masterResultsMatrix_RAS_v2.txt",sep="\t",header=TRUE)

# Constants
types = c("Homdel","Hetloss","Diploid")


shinyServer(function(input, output) {
  
  
  # Show the values using an HTML table
  output$summary <- renderPrint({
    cat(summary(copynumber()),"\n")
    cat(summary(expression()),"\n")
    
  })
  
  
  # Return the requested context
  mutation <- reactive({
    
    if(input$context=="BRAF"){
        mut = braf.mut
    }else if(input$context=="RAS"){
        mut = ras.mut
    }else{
      NULL
    }
    
  })
  
  masterList <- reactive({
    
    if(input$context=="BRAF"){
      master = masterB
    }else if(input$context=="RAS"){
      master = masterR
    }else{
      NULL
    }
    
  })
  
  # Return the requested context
  copynumber <- reactive({
    
    if(input$context=="BRAF"){
      cn = braf.cn2
    }else if(input$context=="RAS"){
      cn = ras.cn2
    }else{
      NULL
    }
    
    # CHEK2
    gene = input$genename
    one.cn = as.numeric(cn[which(rownames(cn)==gene),])
    
  })
  
  # Return the requested context
  expression <- reactive({
    
    if(input$context=="BRAF"){
      exp = braf.exp2
    }else if(input$context=="RAS"){
      exp = ras.exp2
    }else{
      NULL
    }
    
    gene=input$genename
    one.exp = as.numeric(exp[which(rownames(exp)==gene),])
    
  })
  
  methylationRow <- reactive({
    gene=input$genename
    which(rownames(ras.meth2)==gene)
  })
  
  # Return the requested context
  ras.methylation <- reactive({
    
    wh = as.numeric(methylationRow())
    if(length(wh)>0){
      as.numeric(ras.meth2[wh,])
    }else{
      rep(NA,ncol(ras.meth2))
    }
    
  })
  
  braf.methylation <- reactive({
    
    wh = as.numeric(methylationRow())
    if(length(wh)>0){
      as.numeric(braf.meth2[wh,])
    }else{
      rep(NA,ncol(braf.meth2))
    }
    
  })

  samples <- reactive({
    if(input$context=="BRAF"){
      samp = colnames(braf.exp2)
    }else if(input$context=="RAS"){
      samp = colnames(ras.exp2)
    }else{
      NULL
    }
  })
  
  
  
  output$row1plots <- renderPlot({
    
    one.exp = expression()
    one.cn = copynumber()
    mut = mutation()
    gene = input$genename
    expSamps = samples()
    
    numFactors = findInterval(one.cn,c(input$int1[1],input$int1[2]))
    facs = factor(numFactors, labels = types[sort(unique(numFactors))+1])
    
    vec = vector("list",length(unique(facs)))
    vec[[1]] = one.exp[which(facs=="Homdel")]
    vec[[2]] = one.exp[which(facs=="Hetloss")]
    vec[[3]] = one.exp[which(facs=="Diploid")]

    # violin plot
    par(xaxt="n",mfrow=c(1,2))
    vioplot(vec[[1]],vec[[2]],vec[[3]],pchMed=19,colMed="white",rectCol=colors()[222],col=colors()[419],drawRect=FALSE)

    # axis and title and legend
    par(xaxt="t")
    axis(1,at=unique(numFactors)+1,labels=unique(facs))
    title("Violin plot of Expression values")
    
    # now add points: if there is a mutation, that will be plotted separately    
    # mutations
    mm = c()
    if(is.element(gene,mut[,1])){
      ind = which(mut[,1]==gene)
      samps = gsub("-",".",substr(mut[ind,10],1,12))
      mutType = mut[ind,4]
      
      # where to put these on the violin plot
      mm = match(samps,substr(expSamps,1,12))
      mutValues = one.exp[mm]
      xValues = match(facs[mm],types)
      
      mut.ind = as.numeric(mutName[match(mutType,mutName[,1]),2])
      thisMutPch = mutChar[mut.ind[!is.na(mut.ind)]]
      thisMutCol = mutCol[mut.ind[!is.na(mut.ind)]]
    }#end if

    # non-mutated points
    x4points = numFactors
    x4points[mm] = NA
    points(x4points+1,one.exp,col=colors()[130],pch=4)
    
    if(length(mm)>0){
      points(xValues,mutValues,pch=thisMutPch,col=2,bg=colors()[thisMutCol])
      legend("topleft",inset=c(0,0),box.lwd=0,legend=c("missense","nonsense","splice","shift","in_frame"),pch=mutChar,col=2,pt.bg=colors()[mutCol])
    }#end if  
        
    # second plot
    cn.plot = one.cn
    cn.plot[mm] = NA
    plot(cn.plot,one.exp,col=x4points+1,pch=4,ylab="Expression",xlab="Copy number",main="Copy number vs Expression")
    legend("topleft",inset=c(0,0),legend=types,pch=4,col=1:3,,box.lwd = 0)
    if(length(mm)>0){
      points(one.cn[mm],one.exp[mm],pch=thisMutPch,col=2,bg=colors()[thisMutCol])
    }#end if      
  })
  
  output$row2plots <- renderPlot({
    
    one.exp = expression()
    one.cn = copynumber()
    mut = mutation()
    
    # third plot
    par(mfrow=c(1,2))
    plot(density(one.exp,na.rm=TRUE),col=colors()[130],lwd=2,main="Expression density plot")

    # fourth plot
    rasmeth = ras.methylation()
    brafmeth = braf.methylation()
    
    gene=input$genename
    myrow = which(rownames(braf.exp2)==gene)
    brafexp = as.numeric(braf.exp2[myrow,])
    rasexp = as.numeric(ras.exp2[myrow,])
    
    mymin = min(c(brafexp,rasexp)) 
    mymax = max(c(brafexp,rasexp))
    
    myrow2 = which(rownames(braf.cn2)==gene)
    brafcn = as.numeric(braf.cn2[myrow2,])
    rascn = as.numeric(ras.cn2[myrow2,])
    
    numFactors1 = findInterval(brafcn,c(input$int1[1],input$int1[2]))
    numFactors2 = findInterval(rascn,c(input$int1[1],input$int1[2]))
    
    plot(brafmeth,brafexp,ylim=c(mymin,mymax),xlim=c(0,1),col=(numFactors1+1),pch=1,lwd=1,main="Expression vs Methylation",xlab="Methylation",ylab="Expression")
    points(rasmeth,rasexp,col=(numFactors2+1),pch=4,lwd=1)
    legend("topright",legend=c("BRAF+","RAS+"),col=1,pch=c(1,4),box.lwd=0,inset=0)
    legend("topright",legend=types,col=1:3,box.lwd=0,inset=c(0,0.2),lwd=1)
  
  })
  
  sortColumn <- reactive({
    switch(input$scheme,
           "Mann-Whitney test between HETLOSS and DIPLOID expression values"=2,
           "Expression vs copy number distance correlation"=3,
           "Expression unimodality (Hartigan's dip test)"=4,
           "Number of mutations"=5,
           "Expression vs methylation correlation"=6)

  })
  
  sortDirection <- reactive({
    switch(input$scheme,
           "Mann-Whitney test between HETLOSS and DIPLOID expression values"=FALSE,
           "Expression vs copy number distance correlation"=TRUE,
           "Expression unimodality (Hartigan's dip test)"=FALSE,
           "Number of mutations"=TRUE,
           "Expression vs methylation correlation"=FALSE)
    
  })
  
  
  makeSortMas <- reactive({
    master = masterList()
    colNumber = sortColumn()
    direction = sortDirection()
    sortMas = master[order(master[,colNumber],decreasing=direction),]
    rownames(sortMas) = 1:nrow(sortMas)
    
    sortMas[,2] = signif(as.numeric(sortMas[,2]),2)
    sortMas[,3] = round(as.numeric(sortMas[,3]),3)
    sortMas[,4] = round(as.numeric(sortMas[,4]),3)
    sortMas[,6] = round(as.numeric(sortMas[,6]),3)
    
    as.matrix(sortMas)
    
  })
  
  # Show the first "n" genes
  output$view <- renderTable({
    sortMas = makeSortMas()
    head(sortMas, n = input$num2show)
  })
  
  #output$topGene <- renderPrint({
  #  sortMas = makeSortMas()
  #  as.character(sortMas[1,1])
  #})
  
  output$currentGene <- renderText({
    paste(input$genename," in ",input$context,"+ context",sep="")
  })
  
  output$geneStats <- renderTable({
    master = masterList()
    thisrow = which(master[,1] == input$genename)
    thisgene = matrix(NA,nrow=1,ncol=5)
    colnames(thisgene) = c("MW","E.C","Hart","Mut","E.M")
    rownames(thisgene) = input$genename
    thisgene[,1] = signif(master[thisrow,2],2)
    thisgene[,2] = round(master[thisrow,3],5)
    thisgene[,3] = signif(master[thisrow,4],2)
    thisgene[,4] = master[thisrow,5]
    thisgene[,5] = round(master[thisrow,6],3)
    thisgene

  }) 

  
})
