multmerge = function(datalist){
Reduce(function(x,y) {merge(x,y,sort=FALSE)}, datalist)
}