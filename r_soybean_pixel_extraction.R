# Project: 
# Author: Erin Gilbert
# Created:  ~Sept 2017
# Updated: ~ Feb 2019

#Usage:  rscript sds_hyperspec_id.R DARKFILE WHITEFILE RAWFILE
#Purpose: Using thresholds (Mainly NDVI) to create a mask of an image with black being backgrounf (non-plant) and white being plant. These masks are then used to populate a JSON file with the pixel coordinates and cell values of hyperspectral images of plants.


############################################
#####    Define Environment/Packages  ######
############################################

require(sp, lib.loc="/home/hirschc3/gilbe952/source/R/Linux")
require(raster, lib.loc="/home/hirschc3/gilbe952/source/R/Linux")
require(rgdal, lib.loc="/home/hirschc3/gilbe952/source/R/Linux")
require(EBImage, lib.loc="/home/hirschc3/gilbe952/source/R/Linux")
require(rjson, lib.loc="/home/hirschc3/gilbe952/source/R/Linux")
require(foreach, lib.loc="/home/hirschc3/gilbe952/source/R/Linux")

print("libraries loaded")
args <- commandArgs(trailingOnly = TRUE)


############################################
######    Load All the Data    #############
############################################

check_stack_dark<-stack(args[1])
check_stack_white<-stack(args[2])
print("standards loaded")

rawfile<- args[3]
temp2<-strsplit(rawfile, "\\.")[[1]][2]
temp3<-strsplit(temp2, "_")
plot_labels<-strsplit(temp3[[1]][length(temp3[[1]])], "-")[[1]][-1]

hyperspec_stack_raw<-stack(rawfile)
print("raw files loaded")


############################################
######    Perform Calculations    ##########
############################################

wavelengths<-trunc(na.omit(as.numeric(unlist(strsplit(names(hyperspec_stack_raw), split = "X")))))
names(hyperspec_stack_normalized)<-wavelengths


#Use NDvI and PRI to create masks
ndvi_red<-(hyperspec_stack_normalized$X800-hyperspec_stack_normalized$X670)/(hyperspec_stack_normalized$X800+hyperspec_stack_normalized$X670)
pri<-(hyperspec_stack_normalized$X570-hyperspec_stack_normalized$X530)/(hyperspec_stack_normalized$X570+hyperspec_stack_normalized$X530)

ndvi_pri<-ndvi_red
ndvi_pri_label<-thresh(as.array(ndvi_pri), w=300, h=300, offset=0.5)
ndvi_pri_bwlabel<-bwlabel(closing(ndvi_pri_label))
ndvi_pri_bwlabel_table<-sort(table(ndvi_pri_bwlabel), decreasing = TRUE)
final_filter_ndvi_pri<-closing(dilate(rmObjects(ndvi_pri_bwlabel,as.numeric(names(ndvi_pri_bwlabel_table))[(sum(ndvi_pri_bwlabel_table>1000)+1):length(names(ndvi_pri_bwlabel_table))], reenumerate = TRUE)))
#display(final_filter_ndvi_pri)
final_filter_ndvi_pri<-bwlabel(final_filter_ndvi_pri)



# Make sure there are only 6 plants. If more, combine closes objects until 6 reached
center_of_masses_for_plants<-computeFeatures.moment(final_filter_ndvi_pri)[,1:2]
horizontal_midline<-(max(center_of_masses_for_plants[,2])+ min(center_of_masses_for_plants[,2]))/2
print(center_of_masses_for_plants)

#combine_objects
distance_matrix<-as.matrix(dist(center_of_masses_for_plants))

close_objects<-as.matrix(which(distance_matrix< nrow(final_filter_ndvi_pri)/6 & distance_matrix>0, arr.ind = TRUE))[ order(distance_matrix[which(distance_matrix< nrow(final_filter_ndvi_pri)/6 & distance_matrix>0)]),]

for (i in 1:nrow(close_objects))
{
  
  center_of_masses_for_plants<-computeFeatures.moment(final_filter_ndvi_pri)[,1:2]
  if(nrow(center_of_masses_for_plants)==6)
  {
    break
    }else{
    close_objects<-as.matrix(which(distance_matrix< nrow(final_filter_ndvi_pri)/6 & distance_matrix>0, arr.ind = TRUE))[ order(distance_matrix[which(distance_matrix< nrow(final_filter_ndvi_pri)/6 & distance_matrix>0)]),]
    if (close_objects[i,1]!=close_objects[i,2])
    {
      final_filter_ndvi_pri[final_filter_ndvi_pri==max(close_objects[i,])]=min(close_objects[i,])
      close_objects[close_objects==max(close_objects[i,])]<-min(close_objects[i,])
  }else {next}
    }
}


#rename objects
object_names<-c(order(center_of_masses_for_plants[which(center_of_masses_for_plants[,2]<horizontal_midline),1]), order(center_of_masses_for_plants[which(center_of_masses_for_plants[,2]>horizontal_midline),1])+ (nrow(center_of_masses_for_plants)/2))

print(object_names)



############################################
###########    Create Output   #############
############################################



##print the QC images
png(paste("/scratch.",temp2, ".png", sep=""))###file path does start with scratch.global
par(mfrow=c(2,1))
plot(t(ndvi_red), col=gray.colors(15), main=temp2, xaxt='n', yaxt='n', legend=FALSE)
plot(colorLabels(final_filter_ndvi_pri, normalize = T))
dev.off()

# Output a JSON file for each object containing the location and values of the predicted plant pixels

foreach (x=1:length(object_names), .inorder = FALSE) %do% {

  
  result_object<-extract(hyperspec_stack_normalized, which(as.array(final_filter_ndvi_pri)==as.integer(object_names[x])))
  foreach (n=1:ncol(result_object), .inorder=FALSE) %do%{
    result_object[,n]<-(result_object[,n]-cellStats(check_stack_dark[[n]], stat='mean', na.rm=T))/(cellStats(check_stack_white[[n]], stat='mean', na.rm=T)-cellStats(check_stack_dark[[n]], stat='mean', na.rm=T))}
  result_object_json<-NULL
  result_object_json$coordinates<-as.list(data.frame(which(as.array(final_filter_ndvi_pri)==as.integer(object_names[x]), arr.ind = T)))
  result_object_json$hyperspec<-as.list(data.frame(result_object))

  jsonfile<-file(paste(paste("/scratch.",temp2, "Plant", x, sep="-"),"json", sep="."),Sys.getpid(), open = "a" )###file path does not start with scratch.global
  write(toJSON(result_object_json, method="C"), file= jsonfile, append=F)
  close(jsonfile)
}
