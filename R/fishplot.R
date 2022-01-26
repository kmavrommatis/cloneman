
#' infer the clonal model for a cloneobj
#'
#' @param freqs A dataframe with the frequencies of mutations
#' @param founder the clone that is considered to be the founder. If set to NULL the founder will be the clone with the highest CCF or VAF
#' @param model the model of evolution to simulate, (monoclonal|polyclonal)
#' @param type type of frequencies provided (VAF|CCF)

infer.clonal.models.f=function( freqs, founder=NULL, model="monoclonal", type="CCF"){
  # find the clone with the highest vaf
  freqs$cluster = freqs$cluster_name %>% as.factor %>% as.numeric


  if(is.null(founder)){
    freqs[ which(freqs[,2]==max(freqs[,2])), 'cluster'] %>% unique-> founder
    if(length(founder)>1){ founder= min( founder ) }
  }
  res=NA
  if( "VAF"==type){
    res = clonevol::infer.clonal.models(variants=freqs,
                            cluster.col.name="cluster",
                            vaf.col.names =names(freqs)[2:3],
                            subclonal.test="bootstrap",
                            subclonal.test.model="non-parametric",
                            cluster.center="mean",
                            cancer.initiation.model = model,
                            num.boots=1000,
                            founding.cluster=founder,
                            min.cluster.vaf=0.01,
                            sum.p=0.01,
                            alpha=0.1,
                            random.seed=63108) %>%
    clonevol::convert.consensus.tree.clone.to.branch(., branch.scale='none')
  }
  if( "CCF" ==type){
    res = clonevol::infer.clonal.models(variants=freqs,
                                        cluster.col.name="cluster",
                                        ccf.col.names =names(freqs)[2:3],
                                        subclonal.test="bootstrap",
                                        subclonal.test.model="non-parametric",
                                        cluster.center="mean",
                                        cancer.initiation.model = model,
                                        num.boots=1000,
                                        founding.cluster=founder,
                                        min.cluster.vaf=0.01,
                                        sum.p=0.01,
                                        alpha=0.1,
                                        random.seed=63108) %>%
      clonevol::convert.consensus.tree.clone.to.branch(., branch.scale='none')
  }



  mapping=freqs[,c("cluster_name","cluster")] %>% unique
  label=names(freqs)[2]

  res2=lapply( res$models[[label]],function(X){
    Y=merge(X, mapping, by.x="lab",by.y="cluster")
    Y
  })
  res$models[[label]]=res2
  res
}


plot.clonal.models.f=function( res , dir){
  dir.create( dir, recursive=TRUE)
  if(is.null(res)){
    message( "Cannot produce a plot since there are no models")
    return(NULL)
  }
  plot.clonal.models(res,
                     max.num.models.to.plot=1,
                     # box plot parameters
                     box.plot = TRUE,
                     fancy.variant.boxplot.vaf.limits = 100,
                     fancy.boxplot = TRUE,
                     fancy.variant.boxplot.highlight = 'is.driver',
                     fancy.variant.boxplot.highlight.shape = 21,
                     fancy.variant.boxplot.highlight.fill.color = 'red',
                     fancy.variant.boxplot.highlight.color = 'black',
                     fancy.variant.boxplot.highlight.note.col.name = 'gene',
                     fancy.variant.boxplot.highlight.note.color = 'blue',
                     fancy.variant.boxplot.highlight.note.size = 2,
                     fancy.variant.boxplot.jitter.alpha = 1,
                     fancy.variant.boxplot.jitter.center.color = 'grey50',
                     fancy.variant.boxplot.base_size = 12,
                     fancy.variant.boxplot.plot.margin = 1,
                     fancy.variant.boxplot.vaf.suffix = '.VAF',
                     # bell plot parameters
                     clone.shape = 'bell',
                     bell.event = TRUE,
                     bell.event.label.color = 'blue',
                     bell.event.label.angle = 60,
                     clone.time.step.scale = 1,
                     bell.curve.step = 2,
                     # node-based consensus tree parameters
                     merged.tree.plot = FALSE,
                     tree.node.label.split.character = NULL,
                     tree.node.shape = 'circle',
                     tree.node.size = 30,
                     tree.node.text.size = 0.5,
                     merged.tree.node.size.scale = 1.25,
                     merged.tree.node.text.size.scale = 2.5,
                     merged.tree.cell.frac.ci = FALSE,
                     # branch-based consensus tree parameters
                     merged.tree.clone.as.branch = FALSE,
                     mtcab.event.sep.char = ',',
                     mtcab.branch.text.size = 1,
                     mtcab.branch.width = 0.75,
                     mtcab.node.size = 3,
                     mtcab.node.label.size = 1,
                     mtcab.node.text.size = 1.5,
                     # cellular population parameters
                     cell.plot = TRUE,
                     num.cells = 100,
                     cell.border.size = 0.25,
                     cell.border.color = 'black',
                     clone.grouping = 'horizontal',
                     #meta-parameters
                     scale.monoclonal.cell.frac = TRUE,
                     show.score = FALSE,
                     cell.frac.ci = TRUE,
                     disable.cell.frac = FALSE,
                     # output figure parameters
                     out.dir = dir,
                     out.format = 'pdf',
                     overwrite.output = TRUE,
                     width = 4,
                     height = 4,
                     # vector of width scales for each panel from left to right
                     panel.widths = c(3,4,2)
  )
}


#' create a fishplot using the snv events (if available)
#'
#'
#' @param x a cloneobj
#' @param sname name of sample to show on the plot
#' @export
fishPlot=function( x ,sname="sample" ){
  df=getCombinedSNV( x )

  if( max(df$CCF)<0 | max(df$CCF) >1 ){  stop("CCF is not between 0 and 1")}

  label=paste0(X$method_name,".vaf") %>% make.names
  if( "CCF" %in% colnames( S4Vectors::mcols( df ))){
    df$CCF=df$CCF * 100
    ccf=df %>% tibble::as_tibble() %>%
    dplyr::mutate(dummy=CCF) %>%
    dplyr::rename( {{label}} := CCF , cluster_name=clone) %>%
    dplyr::select( cluster_name, !!label, dummy)
    type="CCF"
  }else if("VAF" %in% colnames( S4Vectors::mcols( df ))){
    df$VAF=df$VAF * 100
    ccf=df %>% tibble::as_tibble() %>%
    dplyr::mutate(dummy=VAF) %>%
    dplyr::rename( {{label}} := VAF, cluster_name=clone) %>%
    dplyr::select( cluster_name, !!label, dummy)
    type="VAF"
  }


  y=infer.clonal.models.f(ccf, type=type)
  tryCatch({
    f = clonevol::generateFishplotInputs(results=y)
    fishes = clonevol::createFishPlotObjects(f)
    i=1
    fish = fishplot::layoutClones(fishes[[i]])
    fish = fishplot::setCol(fish,f$clonevol.clone.colors)
    suppressWarnings(
      fishplot::fishPlot(fish,
               shape="spline", title.btm="Patient", cex.title=0.5,
               vlines=1:2, vlab=c(sname), pad.left=0.5)
    )

  },warning=function(w){message(w)}
  , error=function(e){message(e)}

  )
}
