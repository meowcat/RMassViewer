library(RMassBank)
library(mzR)
library(shiny)

# Load data
load("results/cpdList.RData")
files <- list.files(dir.rawfiles, ".mzML", full.names=TRUE) 
load("results/cpdEics.RData")
load("results/cpdMsms.RData")

# Make list of compounds to select from
compoundsSel <- 1:nrow(cpdList)
names(compoundsSel) <- cpdList$CompoundName
# Do the intensity/limit computations like in "plot eics.R"
rtLimits <- do.call(rbind, lapply(cpdEics, function(eic)
{
  do.call(range, lapply(eic, function(e) e$rt ))
}))
intLimits <-do.call(rbind, lapply(cpdEics, function(eic)
{
  do.call(range, lapply(eic, function(e) e$intensity ))
}))
limits <- cbind(rtLimits, intLimits)
# Order intensities
limitOrder <- order(limits[,4], decreasing=TRUE)
# Reorder compound list selection
compoundsSel <- compoundsSel[limitOrder]

filesTitle <- basename(files)

if(!exists("cpdLimits"))
{
  # Make a table for limits loading and saving
  cpdLimits <- cpdList[,c(),drop=FALSE]
  cpdLimits$rtMin <- NA
  cpdLimits$rtMax <- NA
  cpdLimits$rtPeak <- NA
}

server <- function(input, output, session) {
  

  
  # Insert the right number of plot output objects into the web page
  output$plots <- renderUI({
    r.y <- lapply(1:input$panel.y, function(panels.y)
    {
      r.x <- lapply(1:input$panel.x, function(panels.x)
        {
        plotname <- paste("plot",panels.x,panels.y, sep="-")
        clickname <- paste(plotname, "click", sep="-")
        column(width= 12/input$panel.x, plotOutput(plotname, height=input$panel.h, width=input$panel.w,
                                                   clickId=clickname))
      })
      do.call(fixedRow, r.x)
    })
    #do.call(tagList, r.x)
  })
  
  output$plotMapping <- renderUI({
    configFields <- list()
    filesSel <- 1:length(files)
    names(filesSel) <- basename(files)
    filesSel <- c("[none]" = 0, filesSel)
    
    for(panels.y in 1:input$panel.y) {
      for(panels.x in 1:input$panel.x) {
        inputname <- paste("file",panels.x,panels.y, sep="-")
        configFields <- c(configFields,
                          list(selectInput(
                            inputname, inputname, filesSel, input[[inputname]]
                            ))
                          )
      }
    }
    return(do.call(tagList, configFields))
  })


  
  observe({
    input$rtSave
    if(input$rtSave > 0)
      save(cpdLimits, file="results/cpdLimits.RData")
  })

  
  # Update the RT window thing
  observe({
    # Update on changed compound and on loading the RT list
    input$cpd
    input$rtLoaded
    #input$rangeSel
    
    row <- as.integer(input$cpd)
    xlim.eic <- limits[row, c(1,2)]
    #ylim.eic <- limits[row, c(3,4)]
    
   
    if(all(!is.na(cpdLimits[row, c("rtMin", "rtMax", "rtPeak")])))
    {
      print("updating")
      xlim.cpd <- cpdLimits[row, c("rtMin", "rtMax")]
      print("xlim.cpd")
      print(xlim.cpd)
      updateSliderInput(session,"rtLimits",min=xlim.eic[[1]], max=xlim.eic[[2]], value=c(xlim.cpd[1,1], xlim.cpd[1,2]))
      updateSliderInput(session,"rtPeak",min=xlim.eic[[1]], max=xlim.eic[[2]], value=cpdLimits[row, "rtPeak"])  
    }
    else
    {
      print("not updating")
      updateSliderInput(session,"rtLimits",min=xlim.eic[[1]], max=xlim.eic[[2]], value=xlim.eic)
      updateSliderInput(session,"rtPeak",min=xlim.eic[[1]], max=xlim.eic[[2]], value=mean(xlim.eic))  
    }

  })
  

  pool <- reactiveValues()
  pool$rtLoaded <- 0
  
  # load/save compound limits
  observe({
    input$rtLoad
    if(input$rtLoad > 0)
    {
      print("rtload happened")
      file <- "results/cpdLimits.RData"
      cpdLimits <- cpdList[,c(),drop=FALSE]
      cpdLimits$rtMin <- NA
      cpdLimits$rtMax <- NA
      cpdLimits$rtPeak <- NA
      # load if available
      if(file.exists(file))
        load(file)
      cpdLimits <<- cpdLimits
      isolate({pool$rtLoaded <- pool$rtLoaded + 1})
    }
  })
  
  
  
  observe({
    # When limits or peak are changed, set the value in the cpdLimits table
    input$rtLimits
    input$rtPeak
    
    if(all(!is.null(input$rtLimits), !is.null(input$rtPeak)))
    {
      cpd <- as.integer(input$cpd)
      print(cpdLimits[cpd,])
      cpdLimits[cpd, c("rtMin", "rtMax", "rtPeak")] <<-
        c(input$rtLimits, input$rtPeak)
    }
  })
  
  #clickObservers <- reactiveValues()
  clickObservers <- reactiveValues()
  
  # catch the clicks for setting RT range
  observe({
    input$panel.y
    input$panel.x
    # Make an observer for the clickId for every single plot!
    for(panels.y in 1:input$panel.y) {
      for(panels.x in 1:input$panel.x) {
        plotname <- paste("plot",panels.x,panels.y, sep="-")
        clickname <- paste(plotname, "click", sep="-")
        input[[clickname]]
      }
    }
    
    observersList <- list()
    # Make an observer for the clickId for every single plot!
    for(panels.y in 1:input$panel.y) {
      for(panels.x in 1:input$panel.x) {
        plotname <- paste("plot",panels.x,panels.y, sep="-")
        clickname <- paste(plotname, "click", sep="-")
        print(clickname)
        # Observe the clickId for every plot (note: this is a nested observer because otherwise we don't know where the click came from!)
        observersList[[clickname]] <- plotname
      }
    }
    toObserve <- names(observersList)
    observed <- names(clickObservers)
    toKill <- setdiff(observed, toObserve)
    toCreate <- setdiff(toObserve, observed)
    for(clickname in toKill)
    {
      clickObservers[[clickname]]$destroy()
      clickObservers[[clickname]] <- NULL
    }
    for(clickname in toCreate)
    {
      clickObservers[[clickname]] <- observe(substitute({
        input[[clickname]]
        #print(toObserve[[clickname]])
        print(clickname)
        print(input[[clickname]])
      }, list(clickname=clickname)))
    }    
  })
#   observe({clickObservers()})
  
  # The actual plotting function
  observe({
    # replot on change of panel settings
    input$panel.x
    input$panel.y
    for(panels.y in 1:input$panel.y) {
      for(panels.x in 1:input$panel.x) {
        inputname <- paste("file",panels.x,panels.y, sep="-")
        input[[inputname]]
      }
    }
    # replot on selecting a new compound
    input$cpd
    
    # replot when RT limits change
    rangeSel <- input$rangeSel
    rtLimits <- input$rtLimits
    
    if(any(
      is.null(input$panel.x),
      is.null(input$panel.y),
      is.null(input$cpd),
      is.null(input$rangeSel),
      is.null(input$rtLimits)
    ))
      return();
    
    row <- as.integer(input$cpd)
    xlim.eic <- limits[row, c(1,2)]
    ylim.eic <- limits[row, c(3,4)]
    
    if(rangeSel == "window")
      xlim.eic <- rtLimits
    
    for(panels.y in 1:input$panel.y) {
      for(panels.x in 1:input$panel.x) {
        
        local({
          
          inputname <- paste("file",panels.x,panels.y, sep="-")
          plotname <- paste("plot",panels.x,panels.y, sep="-")
 
          if(is.null(input[[inputname]]) || (input[[inputname]] == "0"))
            output[[plotname]] <- renderPlot({
              plot.new()
              title(main=paste(inputname, row))
            })
          else
          {
            file <- as.integer(input[[inputname]])
            output[[plotname]] <- renderPlot({
              par(mar=c(2,1,1,1)+0.1)
              plot.new()
              plot.window(xlim=xlim.eic, ylim=ylim.eic)
              axis(1)
              axis(2)
              eic <- cpdEics[[row]][[file]]
              lines(intensity ~ rt, data=eic, col=fileColors[[file]])
              msms <- cpdMsms[[row]][[file]]
              for(specs in msms)
              {
                if(specs$foundOK == TRUE) {
                  segments(specs$childHeaders$retentionTime, ylim.eic[[2]], specs$childHeaders$retentionTime, ylim.eic[[2]]/2, col="red")
                  #         lines(basePeakIntensity ~retentionTime, data = specs$childHeaders, 
                  #               col="red", type='h' )
                } # endif(specs$foundOK)
              } # endfor(specs in msms) 
              legend("topleft", filesTitle[[file]],  bty="n")
              
            }) # renderPlot
          } # else[if file != 0]
        
        }) #local
      } # for panels.x
    } # for panels.y

  }) # observe  


  # save and load settings
  observe({
    input$settingsSave
    
    if(input$settingsSave > 0)
    {
      isolate({
        if(file.exists("results/RmvSettings.RData"))
          load("results/RmvSettings.RData")
        else
          RmvSettings <- list()
        newSettingsName <- input$settingsName
        newSettings <- list(
          panel.x = input$panel.x,
          panel.y = input$panel.y,
          panel.w = input$panel.w,
          panel.h = input$panel.h,
          rangeSel = input$rangeSel
        )
        for(panels.y in 1:input$panel.y) {
          for(panels.x in 1:input$panel.x) {
            inputname <- paste("file",panels.x,panels.y, sep="-")
            newSettings[[inputname]] <- input[[inputname]]
          }
        }
        RmvSettings[[newSettingsName]] <- newSettings
        updateSelectInput(session,"settingsList", choices=names(RmvSettings), selected=newSettingsName)
        save(RmvSettings, file="results/RmvSettings.RData")
      })
    }
  })
  # save and load settings
  observe({
    input$settingsLoad
    
    if(input$settingsLoad > 0)
    {
      isolate({
        load("results/RmvSettings.RData")
        newSettingsName <- input$settingsList
        settings <- RmvSettings[[newSettingsName]]
        updateSliderInput(session, "panel.x", value=settings$panel.x)
        updateSliderInput(session, "panel.y", value=settings$panel.y)
        updateSliderInput(session, "panel.w", value=settings$panel.w)
        updateSliderInput(session, "panel.h", value=settings$panel.h)
        updateRadioButtons(session, "rangeSel", selected=settings$rangeSel)
        
        for(panels.y in 1:settings$panel.y) {
          for(panels.x in 1:settings$panel.x) {
            inputname <- paste("file",panels.x,panels.y, sep="-")
            updateSelectInput(session, inputname, selected=settings[[inputname]])
          }
        }
      })
    }
  })

# delete settings
observe({
  input$settingsDel
  
  if(input$settingsDel > 0)
  {
    isolate({
      load("results/RmvSettings.RData")
      newSettingsName <- input$settingsList
      RmvSettings[[newSettingsName]] <- NULL
      updateSelectInput(session,"settingsList", choices=names(RmvSettings))
      save(RmvSettings, file="results/RmvSettings.RData")
      
    })
  }
})
  
  if(file.exists("results/RmvSettings.RData"))
  {
    load("results/RmvSettings.RData")
    updateSelectInput(session,"settingsList", choices=names(RmvSettings))
  }
}

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      selectInput("cpd", "Compound", compoundsSel, compoundsSel[[1]]),
      radioButtons("rangeSel", "Range",
                   c(               "full range" = "full",         "RT window" = "window"),"full"),
      conditionalPanel(condition = "input.rangeSel == 'window'",
                       sliderInput("rtLimits", "RT range", 0,999,c(0,999))),
      sliderInput("rtPeak", "Peak RT", 0,999,500),
      actionButton("rtSave", "save RT limits"),
      actionButton("rtLoad", "load RT limits")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("View",
                 uiOutput("plots")
                 ),
        tabPanel("Setup",
                 actionButton("settingsLoad", "load layout"),
                 actionButton("settingsDel", "delete layout"),
                 selectInput("settingsList","Layout", c()),
                 actionButton("settingsSave", "save layout"),
                 textInput("settingsName", "Name",""),
                 sliderInput("panel.x", "Horiz. panels", 1,8,3),
                 sliderInput("panel.y", "Vert. panels", 1,8,3),
                 sliderInput("panel.w", "Panel width", 250,600,300),
                 sliderInput("panel.h", "Panel height", 250,600,300),
                 uiOutput("plotMapping")
                 )
        )
      )
  )
)

shinyApp(ui = ui, server = server)


