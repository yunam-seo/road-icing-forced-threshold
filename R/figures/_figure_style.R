# Opens the graphics device for the figure scripts (PNG by default, PDF when FIGURE_VECTOR=1) and adds panel letters.
# Sourced by the figure scripts 62-67.

.pl_env<-new.env(); .pl_env$q<-character(0)
panel_letters<-function(l){ .pl_env$q<-l; invisible(NULL) }
setHook("plot.new",function(){
  q<-.pl_env$q
  if(length(q)){
    text(grconvertX(0.012,"nfc","user"),grconvertY(0.992,"nfc","user"),q[1],font=2,cex=1.35,adj=c(0,1),xpd=NA)
    .pl_env$q<-q[-1] }
},action="replace")
figure_device<-function(path,width,height,res,pointsize){
  if(!identical(Sys.getenv("FIGURE_VECTOR"),"1")){
    png(path,width=width,height=height,res=res,pointsize=pointsize)
    return(invisible(path))
  }
  vd<-file.path(dirname(path),"VECTOR"); dir.create(vd,showWarnings=FALSE,recursive=TRUE)
  p<-file.path(vd,sub("[.]png$",".pdf",basename(path)))
  ok<-tryCatch({grDevices::cairo_pdf(p,width=width/res,height=height/res,
                                     pointsize=pointsize,onefile=FALSE); TRUE},
               error=function(e) FALSE)
  if(!ok) pdf(p,width=width/res,height=height/res,pointsize=pointsize,onefile=FALSE)
  message("  [vector] ",p)
  invisible(p)
}
