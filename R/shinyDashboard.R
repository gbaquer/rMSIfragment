#Shiny dashboard
is<-(1:length(pks$mass))
names(is)<-as.character(round(pks$mass,2))
ui <- dashboardPage(
  dashboardHeader(title = "rMSIfragment explorer"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("MSI Dataset", tabName = "dataset", icon = icon("shapes")),
      menuItem("Annotations", tabName = "annotations", icon = icon("th-list")),
      menuItem("Explorer", tabName = "explorer", icon = icon("project-diagram"))
    )
  ),
  dashboardBody(
    tabItems(
      # Sample Content
      tabItem(tabName = "dataset",
              fluidRow(plotOutput("meanSpectra")),
              fluidRow(
                column(4,selectInput("i1","Ion 1",is),plotOutput("ion1")),
                column(4,selectInput("i2","Ion 2",is),plotOutput("ion2")),
                column(4,selectInput("i3","Ion 3",is),plotOutput("ion3"))
              )
      ),
      #Annotation Content
      tabItem(tabName = "annotations",
              fluidRow(
                column(12,div(DTOutput("res"), style = "overflow-y: auto;"))
              ),
              downloadButton("downloadData", "Download .csv")
      ),

      # Explorer Content
      tabItem(tabName = "explorer",
              fluidRow(column(2,selectInput("mz","Choose mz",is)),column(2,selectInput("n","Page",1:5)),column(3,radioButtons("nID","Choose network",choices=1:5,inline=T))),
              fluidRow(div(style = 'overflow-y: scroll',plotOutput("networks"))),
              fluidRow(plotOutput("plots"))
      )
    )
  )
)

server <- function(input, output) {

  output$res <- renderDT({datatable(simply(res),filter = "top")%>%
      formatRound(columns = c(1,2,8), digits = 2)})

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("annotation-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(res, file)
    }
  )
  output$ion1 <- renderPlot({
    ggplot_peak_image(one_pks(pks,4),as.numeric(input$i1))

  })
  output$ion2 <- renderPlot({
    ggplot_peak_image(one_pks(pks,4),as.numeric(input$i2))
  })
  output$ion3 <- renderPlot({
    ggplot_peak_image(one_pks(pks,4),as.numeric(input$i3))
  })

  output$meanSpectra <-renderPlot({
    ggplot_mean_spectra(one_pks(pks,4),input$i1,input$i2,input$i3)
  })

  k<-reactive({plot_networks(one_pks(pks,4),res5ppm,pks$mass[as.numeric(input$mz)])})

  output$networks <- renderPlot({
    display_networks(k(),as.numeric(input$n),5)
  })

  output$plots <- renderPlot({
    plot_ions(one_pks(pks,4),res5ppm,pks$mass[as.numeric(input$mz)],as.numeric(input$nID)+(as.numeric(input$n)-1)*5)
  })
}
plot_networks<-function(pks,res,mz,topN=NA,score="lipid_occurences*(1+correlation)"){
  res<-res[order(with(res,eval(parse(text=score))),decreasing=T),]
  i<-which.min(abs(mz-pks$mass))
  x<-subset(res,mz_i==i)
  if(!is.na(topN)&nrow(x)>0){
    x<-x[1:min(nrow(x),topN),]
  }
  k<-list()
  for(i in 1:nrow(x)){
    y<-subset(res,lipid_id==x$lipid_id[i])
    z<-matrix(0,nrow(y)+1,nrow(y)+1)
    n=with(rbind(x[i,],y),paste(round(experimental_mz,2),"\n[",adduct," ",fragmentation,"]",sep=""))
    n[1]<- paste(x$abbreviation[i],x$formula[i],n[1],paste("LO:",x$lipid_occurences[i],"C:",round(x$correlation[i],2)),sep="\n")
    z[1,]=1
    net = network(z, directed = T)
    k[[i]]<-ggnet2(net,mode="kamadakawai",label=n,shape=15,label.size=3.5,size=20,color=1+c(0,as.numeric(y$fragmentation=="")))+
      scale_color_brewer(palette="Set3",labels = c("Annotation for current mz","Fragment","Parental ion"))+theme(aspect.ratio = 1)
  }
  return(k)
}
display_networks<-function(k,i,n=5){
  p<-grid.arrange(grobs=k[1:length(k)%in%(((i-1)*5+1):(i*n))],nrow=1)
  return(p)
}
plot_ions<-function(pks,res,mz,j,topN=NA,score="lipid_occurences*(1+correlation)"){

  res<-res[order(with(res,eval(parse(text=score))),decreasing=T),]
  i<-which.min(abs(mz-pks$mass))
  x<-subset(res,mz_i==i)
  if(!is.na(topN)&nrow(x)>0){
    x<-x[1:min(nrow(x),topN),]
  }
  y<-rbind(x[j,],subset(res,lipid_id==x$lipid_id[j]))
  p<-grid.arrange(grobs=lapply(1:nrow(y),function(i)ggplot_peak_image(pks,y$mz_i[i])))
  return(p)
}
#Uncomment to run
#shinyApp(ui, server)
