########### TABLE OF CONTENTS ############################
#0.  Define all libraries and functions
#1.  File Handling - manually enter file name as variable f 
#           (USER INPUT: inpu file name or construct loop)
#2.  Extract meta-tags from thermal vid file
#3.  File information for binary read in  
#           (USER INPUT: check width and length)
#4.  Read first 10,000,000 bytes Binary data input 
#5.  Find everywhere in vector where resolution info is: i.e. 320x240
#6a. Read header information in with readBin and seek 
#6b. Read/seek and fill the frame with data.
#7.  Extract time points from the raw headers
#8.  Convert raw data in alldata variable to temperature
### -- Below here contents should be the same as found in Flir files ---
#9.  Image variability analysis to estimate activity difference frame by frame.
#10. Convert activity data into cumulative sums 
#11. Calculate slope of the cumulative function every N samples 
#          (USER INPUT: # of Samples to avg)
#12. Put final activity data into .csv files 
#         (Double check folder location that was set in section 1)
#13. Plot the thermal images.


#0. Define Libraries and functions ####
library("lattice")
library("fields")
library("animation")
library("Thermimage")
library("raster")

#----- Set the colour palette to be used 
# for use with image.plot and other rasterised thermal images
thermal.palette<-palette.choose("ironbow")  # can choose form "flir","ironbow"...need to add others


#1. File Handling - manually enter file name as variable f ####
# define the file output names and locations

f<-NULL

mainDir<-"/Users/GlennTattersall/Desktop/thermalvids/"
cat("Files will be loaded from: ", mainDir)
cat("\n")

outputDir<-"/Users/GlennTattersall/Desktop/output/"
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
cat("Output files will be saved in: ", outputDir)
cat("\n")

vidDir<-paste(outputDir,"vid_output",sep="")
dir.create(vidDir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
cat("Video output files will be saved in: ", vidDir)
cat("\n")

setwd(mainDir)       
getwd() 
#show working directory


l.files<-list.files(mainDir,all.files=FALSE,full.names=FALSE,recursive=TRUE,
                    include.dirs=TRUE, pattern = "\\.rtv$")

fn<-l.files[1]   # rem this out and reconstruct the loop below for analysing folders of files

#for(fn in l.files[1:49])
#{
  
cat("Analysing file:" ,fn, sep=" ")
cat("\n")
cat("File #:", which(fn==l.files), "of", length(l.files), sep=" ")
cat("\n")


f<-paste(mainDir,fn,sep="")
  
fname<-unlist(strsplit(fn,"[/]"))
fname<-fname[length(fname)]
# Remove the 'directory information' ahead of the actual filename
# fname is used for constructing the f.root used in output files

getwd() #show working directory
setwd(mainDir)              
  
f.root<-substr(fname,1,nchar(fname)-4)     
# text output file name root, without .rtv/seq
  
outputID_Dir<-paste(vidDir,"/", f.root, sep="")
dir.create(outputID_Dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
cat("This file's video output will be saved in: ")
cat(outputID_Dir)
cat("\n")


txt.outDir<-paste(outputDir,"txt_output",sep="")
# data arising from this script will be saved here
dir.create(txt.outDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
cat("This file's text output files will be saved in: ")
cat(txt.outDir)
cat("\n")


f.raw.activity<-paste(txt.outDir, "/", f.root,"_raw.activity.csv", sep="")
f.resamp.activity<-paste(txt.outDir, "/", f.root,"_resamp.activity.csv", sep="")
  

#2. Extract meta-tags from thermal vid file ####
# this is a placeholder section now.  At the moment, there is no meta-tag data from Mikron .rtv files 
# since the file structure is very simple and contains no information on the camera.
# Perhaps in the future there will be.


#3. File info for read-in ####
finfo = file.info(f)
finfo$size  # how many bytes in the file
byte.length<-2  # how many bytes make up an element
no.elements<-finfo$size/byte.length 
# this should tell you how many total elements make up the inputted data
op<-options(digits.secs=3)
save.time<-strptime(finfo$mtime,"%Y-%m-%d %H:%M:%OS")
# time at which the file was saved.  this should correspond to the last frame of the file
# this does the same: save.time<-finfo$mtime
# DO NOT TRUST: This is a "file modified time" value and may not correspond to the time on that the
# last frame was recorded
w<-320    # set image width for MikrospecRT image
h<-240    # set image height for MikrospecRT image
l<-w*h



#4.  Read initial 10,000,000 bytes binary Data in ####
# Open file for populating alldata with imaging data
dur<-20   # duration of frame header is 20 (actually 21, but this is added to start element)
s<-1
e<-s+dur
footerlen<-810 # last 810 bytes are also a footer in rtv file
headerlen<-42
# read in file, with 2 bytes per element
to.read = file(f, "rb") # set to.read file.  rb means read binary
alldata<-readBin(to.read, "integer", n=10000000, 
                 size=byte.length, endian = "little", signed=FALSE)
#no.con<-round(length(alldata)/l) # number of contaminating elements (i.e. headers/frames) in image
close(to.read)


#5. Find everywhere in vector where resolution info is: i.e. 320x240 ####
# this allows to find the start of each header and data frame 
# note: user must set l and w at the top of this script          
fid<-c(2,h) # if using little endian
# trouble with very large files, so try a small portion and predict the remaining 
# wh.locates

if(length(alldata)>=1000000)
{
  system.time(wh.locate<-locate.fid(fid,alldata,long=TRUE))
  # try wh.locate on a small chunk of all data
  diff.wh.locate<-diff(wh.locate) 
  # difference calc should yield a repeating pattern of header then frame
  gaps<-unique(diff.wh.locate)
  # if the pattern is simple, this value should be 2
  no.unique.locates<-length(gaps) 
  # reconstruct wh.locate from scrap, starting only with wh.locate[1]
  repeats<-2*trunc(finfo$size/2/(l+wh.locate[1]+gaps[1]))
  # how many repeats required to create fill whole file
  wh.locate<-cumsum(as.numeric((c(wh.locate[1],rep(gaps,repeats)))))
  # cumulative sum up the 1st locate and repeate gaps
  
  wh.locate<-wh.locate[-c(which(finfo$size/2-wh.locate<l))]
  # remove any header start that is too close to the end of the file to be a potential frame
  
  if (length(which(wh.locate>(finfo$size/byte.length)))>0)
    {
      wh.locate<-(wh.locate[-c(which(wh.locate>(finfo$size/byte.length)))])
      # remove any locates that go beyond the length of the data file after import    
    }
  
}

if(length(alldata)<1000000)
{
  # much faster without calling my function:
  system.time(fid1.locate<-which(alldata==fid[1]))
  fid1.locate.adjacent<-fid1.locate+1
  fid2.locate<-which(alldata[c(fid1.locate.adjacent)]==fid[2])
  wh.locate<-fid1.locate[fid2.locate]
  
  diff.wh.locate<-diff(wh.locate) 
  # difference calc should yield a repeating pattern of header then frame
  gaps<-unique(diff.wh.locate)
  # if the pattern is simple, this value should be 2
  no.unique.locates<-length(unique(diff.wh.locate)) 
  # if the pattern is simple, this value should be 2
}

header.l<-diff.wh.locate[length(diff.wh.locate)]
# the approximate length of the header, as defined by the last element in the diff.wh.locate vector
# the frame data stream commences 9 elements after the resolution info
res2fram<-9
# Below define the start indices for the headers (h.start) and frames (f.start)
# check if the first location of the resolution info is a small value or not
if(wh.locate[1]<header.l & no.unique.locates==1)  # .rtv files appear to be formatted this way  
{
  h.start<-wh.locate-12
  h.end<-wh.locate+res2fram-1
  f.start<-wh.locate+res2fram
  f.end<-f.start+l-1
}
# each f.start should correspond to the start of a frame in the video file
# from my impression, the location of the header is 9 elements in front of the first pixel value
# above, res2fram is set to be 9



# 6a. Read header information in with readBin and seek ####
# since we have defined the h.starts and h.ends above, in principle this can be used
# to open the "to.read" variable at only those locations
# For raw data imported 1 byte at a time, it is possible to extract frame times that are held in every header
# In fact, the first 22 bytes of an .rtv file is the data and time of the first frame.
# 153619 bytes later, the next time stamp can be found Empirically determined # of bytes where repeating
# but can be calculated as well. Each header is 42 bytes long.  
# The date occupies the first 22 bytes of the header.  To skip to the next header, then move:
# 42+2*(320*240) = 153642 bytes
# To estimate the # of frames in a file, take (filesize-footerbytes)/(2*320*240+42)
# Round the number to ensure an integer value  
# footerlen and headerlen are defined earlier in file: 810 and 42 bytes each
timeRaw<-NULL
dateRaw<-NULL
to.read <- file(f, "rb") # set to.read file.  rb means read binary
for(i in 1:length(h.start))
{
  t.s<-(i-1)*153642      # start bytes for where date/time stamp is
  seek(to.read,where=t.s,origin="start")
  timeRaw<-readBin(to.read, raw(), n=22, size=1, endian = "little", signed=FALSE)
  dateRaw<-cbind(dateRaw,timeRaw)
  cat('\r', paste("Importing Time stamps from header # ", i, " of ", length(h.start), sep=""))
  Sys.sleep(0.0005)
}
close(to.read)
cat('\n')
colnames(dateRaw)<-NULL


# 6b. Read/seek and fill the frame with data.   ####
# using a data.frame filling method.  should be more healthy for RAM and recalling
# frame by frame data

fram<-NULL
flagged.fram<-NULL
alldata<-NULL
to.read <- file(f, "rb")
system.time(for (i in 1:(length(f.start)))
{
  seek(to.read,where=(f.start[i])*byte.length-2,origin="start")
  fram<-readBin(to.read, integer(), n=l, 
                size=byte.length, endian = "little", signed=FALSE)
  #colnames(fram)<-"Frame"
  if(length(which(fram==0))>0) flagged.fram<-c(flagged.fram,i)
  fram[which(fram==0)]<-1
  if (i==1) alldata<-data.frame(fram,check.names=TRUE)
  if(i>1) alldata<-data.frame(alldata,fram,check.names=TRUE)
  cat('\r', paste("Importing Frame # ", i, sep=""))
}
)
close(to.read)
cat('\n')
cat(paste("The following frames were flagged for errors:"))
cat(flagged.fram)
flagged.fram.index<-rep(0,length(f.start))
flagged.fram.index[flagged.fram]<-1
cat('\n')


#system.time(h.ind<-unlist(Map(':',h.start, h.end)))
#system.time(f.ind<-unlist(Map(':',f.start, f.end)))


#7. Extract time points from raw inputs ####
dateChar<-apply(dateRaw,2,rawToChar)
extract.times<-strptime(dateChar,format="%m/%d/%y %H:%M:%OS")
# time values expressed in posixct format
dateOriginal<-extract.times[1]

# Extract the time values from the headers:
diff.times<-as.numeric(diff(extract.times))
cum.times<-cumsum(diff.times)
f.times<-0
f.times<-c(f.times,cum.times)
options(digits.secs=3)
est.times<-extract.times[length(extract.times)]-rev(f.times)-f.times[1]
# real times in Date:Time format based on finfo$mtime which should correspond to the time
# the file was saved (i.e. save.time<-finfo$mtime.  Only problem with this is that if you 
# edit or re-save a .seq file, then this will not work.  
# So the est.times will essentially be the same as extract.times, but I've included it here
# for posterity and comparison to other scripts


#8. Convert raw data held in alldata variable to temperature ####
# First estimate a true 'raw' value from Mikron's raw, since Mikron's raw is simply 10*Temperature (degK)
# no need for special Mikron calibration constants since the raw.est will be converted into an adjusted
# temperature measurement. However, must multiply the alldata by 0.1 and then substract 273.15 in order to express the Mikron 'raw'
# data in appropriate temperature units.
# raw.est is then our best estimate of the raw, non-linearised units of the detector, which is done by 
# setting all E,T and tau values to unity.
#system.time(rangemikron<-range(alldata))
mikronlookup<-(1:7730)
# a false mikron vector starting at 1 (a value that mikron does not produce) and goes to 7730
# 7730 corresponds to an estimate Mikron temperature of 500C --> 773K --> 7730 raw
rawlookup<-temp2raw((0.1*mikronlookup-273.15),1,0,20,20,20,1,50)
# set rawlookup to values that would resemble those at calibration

#raw.est<-temp2raw(alldata*0.1-273.15,1,0,20,20,20,1,50)

# once you have raw.estimate, then use parameters relevant to experiment to estimate temperature.
# The same calibration constants (defaults are fine) used to generate the raw.est should be used to generate
# the temperature estimates, but technically this is not being 'calibrated', but using the functions to
# tranform data according to the physical conditions of the actual measurements.
# For our lab measurements, here are typical conditions:
# Distance 0.5 m, Feather/body emissivity 0.96
# Optical window transmittance 0.96, temperatures set to 20C, RH set to 50%
# Technically the RTemp, ATemp and IRWTemp should be the incubator temperature.
# Since Mikron's software does not account for atmospheric 
templookup<-raw2temp(rawlookup,0.96,0.4,20,20,20,0.96,50)

    # to use templookup-->  templookup[alldata[1:10,1]]

#alltemperature<-raw2temp(raw.est,0.96,1,20,20,20,0.96,50)


#alltemperature<-raw2temp(alldata,0.96,100,20,20,20,0.96,50, PR1=21106.77,PB=9758.743281,
#                       PF=29.37648768,PO=1278.90,PR2=0.03766376)
# Note: since Mikron stores 'raw' data as T, of a blackbody in deg K, the calibration I have above
# is to be compatible with FLIR's algorithm.  This means that the raw 'alldata' signal must be 
# multiplied first by 0.1 in order to convert the raw signal into bb, Temperature

rangetemp<-NULL
rangetemp[1:length(f.start)*2]<-NA
j<-seq(1,length(f.start)*2,2)
system.time(
  for (i in 1:length(f.start))  
  {
    rng.temp<-range(templookup[alldata[1:(l-1),i]], na.rm=TRUE)
    rangetemp[j[i]:(j[i]+1)]<-rng.temp
    cat('\r',paste("Calculating Frame #", i, "   Temperatures-> Min:", round(rng.temp[1],digits=1), "Max:", 
                   round(rng.temp[2], digits=1), "and Mid:", round(sum(rng.temp)/2,digits=1), "     ", sep=" "))
  }
) 

mintemp<-rangetemp[c(seq(1,length(rangetemp),2))]
maxtemp<-rangetemp[c(seq(2,length(rangetemp),2))]
midtemp<-(maxtemp+mintemp)/2
allmintemp<-median(mintemp, na.rm=TRUE)
allmaxtemp<-median(maxtemp, na.rm=TRUE)
allmidtemp<-median(midtemp, na.rm=TRUE)
cat('\n')


#9. Image variability analysis to estimate activity difference frame by frame. ####
avg<-NULL
sdev<-NULL
d<-rep(0,l)
for (i in 2:length(f.start)) {       
  # use this for loop to open file by file and perform diff analysis
  prev.i<-i-1             # index for the previous 'frame'
  #pfs<-i                  # previous frame start point
  #pfe<-pfs+l-1            # previous frame end point
  #cfs<-i                  # current frame start point
  #cfe<-cfs+l-1            # current fram end point
  dd<-abs(alldata[,i]-alldata[,prev.i])
  # calculate the diff between individual elements (i.e. pixels) between frames
  # store in temporary variable dd
  if(abs(mean(dd)-mean(d))<5000) d<-dd
  # check difference in the frame differences to make sure it is not a calibration sweep
  # if it is a small change (i.e. less than 5000), then set d equal to dd.
  # otherwise d will remain the same as the previous d value
  # print(c(i))
  cat('\r',paste("Calculating differences between Frame # ",i-1, " and ", i,
                 which(d==NA), sep=""))
  #Sys.sleep(0.01)
  #if(i>1 & sd(d)==0)  
  # break out of for loop the first instance of a frame freeze occurs (common in Mikron
  # when there is a memory problem during data recording and all images are identical)
  #{
  #  break
  #}
  av<-mean(d,na.rm=TRUE)           # average of the |inter-frame difference|
  sde<-sd(d,na.rm=TRUE)            # sd of the |inter-frame difference|
  avg<-c(avg,av)        # bind the current av to previous avg to create a vector of avg abs(diff) values
  sdev<-c(sdev,sde)     # bind the current sde to previous sdev to create a vector of sd of abs(diff) values
}
cat('\n')



#10. Convert activity data into cumulative sums  ####
# The avg and sdev values are simply frame by frame assessments of change. If the frame rate is fast,
# then changes between frames will be small. To turn this into more meaningful data involves calculating a 
# cumulative sum function of each of these vectors. 
# This is referred to as a cumulative summation function, where the 1st value
# is unchanged, but each successive element in the vector is the sum of the previous 
# and itself. Total length of this vector is unchanged, since the first element is simply itself.
# The second element in the cumsum function is the sum of the second + first element of the input.
# The third element in the cumsum function is the sum of the third + second element of the input.
# and so on...
activity<-NULL
Fs<-1/mean(diff(f.times)) # sample rate in frames per sec (images captures every 1/Fs seconds)
activity<-cbind(avg,sdev)
no.cols<-ncol(activity)
colnames(activity)<-c("Avg","SD")
rownames(activity)<-c(seq(1,nrow(activity),1))
cumsum.act<-apply(activity,2,cumsum)  
# calculate the cumulative sum for each of the columns in matrix 'activity'
cumsum.act<-cumsum.act-cumsum.act[1,] 
# subtract element 1 from each. Tidies things up so that cumum.act starts at a value of zero
df.time<-seq(1/Fs,nrow(cumsum.act)/Fs,1/Fs)
# construct proper time index to accompany the cumulative sum data
df.times<-est.times[2:length(est.times)]
# construct actual time index from the est.times variable used earlier during header extraction
plot.new()
matplot(df.time, cumsum.act, type = "l",  lty=c(1,1), lwd=1, xlab="Time (s)", ylab="Cumulative Sum",
        ylim=c(min(cumsum.act), max(cumsum.act)), xlim=c(min(df.time),max(df.time)))
legend("topleft", bty = "n", c("Avg", "SD"), lty=c(1,1), col=c("black", "red"))
# plot the cumulative function 


#11. Calculate slope of the cumulative function every N samples ####
samples<-4 # set this value to the number of sample points over which to calculate avg(slope) data 
# if your sample rate is 1 frame per second, then this is simply a value every 20 seconds
# if your sample rate is 33 frames per second, then this is a value every 20/33 = 0.6 seconds
slp.av<-slopeEveryN(cumsum.act[,1],samples)
slp.sd<-slopeEveryN(cumsum.act[,2],samples)
slp.sd<-slp.sd[,-1]                  # remove the 1st column from slp.sd
slp<-cbind(slp.av,slp.sd)            # bind slp.av and slp.sd: slopes taken every n samples
#slp[,2]<-slp[,2]-slp[1,2]            # subtract first element in order to have slopes start ~0
#slp[,3]<-slp[,3]-slp[1,3]            # subtract first element in order to have slopes start ~0
len<-round(length(df.times)/samples)  # how many 'rows' make up the resampled data
#time.s<-seq(from = samples/Fs, to = samples/Fs*len, by = samples/Fs)
#time.m<-time.s/60                     # convert times to minutes
# df time variables refers to the frame differences. 
# time variables refer to absolute clock time
# time.s or time.m refer to elapsed time in seconds or minute
dfsec<-seq(from = 1/Fs, to = length(df.times)/Fs, by = 1/Fs)
dfmin<-dfsec/60                     # convert times to minutes
samp.times<-meanEveryN(df.times,samples)
#samp.times<-as.numeric(samp.times[,-1])
samp.times<-as.POSIXct(samp.times, origin="1970-01-01")  # extracted real times at sample rate
#samp.sec<-meanEveryN(dftime.s,samples)
samp.sec<-meanEveryN(dfsec,samples)
#samp.sec<-(samp.sec[,-1])
samp.min<-samp.sec/60
matplot(samp.min, slp[,2:3], type = "l",  lty=c(1,1), lwd=1, xlab="Time (min)", ylab="Slope of cumulative sum",
        ylim=c(min(slp[,2:3]), max(slp[,2:3])), xlim=c(min(dfmin),max(dfmin)))
legend("topleft", bty = "n", c("Avg", "SD"), lty=c(1,1), col=c("black", "red"))


#12. Put final activity data into csv files ####
nf<-length(f.start) # quick count of the number of frames to be exported 
# although this value (nf) is the total frame #.  for the analysis of frame difference analysis: nf-1
activity.output<-data.frame(2:nf, flagged.fram.index[2:nf], extract.times[2:nf], est.times[2:nf], f.times[2:nf], 
                            activity, cumsum.act, mintemp[2:nf],
                            maxtemp[2:nf], midtemp[2:nf])
colnames(activity.output)<-c("Frame #", "Flagged Frame", "Ext. Times", "Est. Times", "Elapsed Time(s)", 
                             "Average", "SD", "CumulSumAvg", "CumulSumSD", "Min Temp", 
                             "Max Temp", "Mid Temp")
resam.output<-data.frame(samp.times,samp.min,samp.sec,slp)
colnames(resam.output)<-c("Samp Times", "Elapsed Time (min)", "Elapsed Time(s)","Sample #", "Slope Avg","Slope SD")
write.csv(activity.output, f.raw.activity)
write.csv(resam.output, f.resamp.activity)


#13. Plot the thermal images. ####
# This requires the "animation" package, which is not available in ready form on CRAN. 
# Therefore, you must download the .tar file to your computer and install from there.
# install.packages("~/Downloads/animation_2.3.tar", repos= NULL, type="mac.binary.mavericks")  
# library("animation") need to download seewave tar from website:
# http://www.renevolution.com/how-to-install-ffmpeg-on-mac-os-x/
# Save image sequence as video, with image number of and time. saveVideo
#Or use saveHTML function
#Then Open images in R and take thermometrics. (not implemented yet)
plot.new()
des = c(f.root)
setwd(outputID_Dir)
subtitle<-paste(as.character(extract.times)," Frame #: ")
minrangeset<-21
minrangeset<-41
#minrangeset<-0.1*rangemikron[1]-273.15
#maxrangeset<-0.1*rangemikron[2]-273.15
system.time(
  saveHTML(  
  #saveGIF(  
  #saveVideo( 
    for(i in 1:length(f.start))
    {
      ani.options(interval = 0.08, ani.width=640)
      #ani.options(interval = 1/Fs)
      sub.title<-paste(subtitle[i],i)
      d<-matrix(rev(templookup[alldata[,i]]), nrow=w, ncol=h)
      d[which(d>maxrangeset)]<-maxrangeset
      d[which(d<minrangeset)]<-minrangeset
      # use the templookup variable, to reduce the requirement for alltemperature
      # populate temperature matrix in reverse in order to flip image before plotting 
      d<-flip.matrix(d)
      #dr<-raster(rotate90.matrix(d), xmn=0, xmx=640, ymn=0, ymx=480) # convert matrix to raster
      #  ds<-(d-mintemp)/(maxtemp-mintemp) 
      #par(pin=c(4,3)) 
      #plot(dr,add=FALSE, useRaster=TRUE, box=FALSE, axes=FALSE, main=paste(maintitle[i],i), 
      #     xlim=c(0,640), ylim=c(0,480), zlim=c(mintemp,maxtemp), 
      #     legend.shrink=0.8, legend.width=1, col=thermal.palette)
      par(pin=c(6,4.5),cex.sub=1.5, cex.main=1.5)
      image.plot(d, useRaster=TRUE, bty="n", col=thermal.palette, main=des,
                 xlab="", ylab="", xaxt="n", yaxt="n", 
                 zlim=c(minrangeset,maxrangeset),
                 legend.shrink=0.85, legend.cex=0.85, asp=h/w)
      title(sub=sub.title, adj=0.5)
      title(cex=1, adj=0.5, sub=paste("Min: ", format(round(mintemp[i], 1),nsmal=1), "C ",
                                      " Max: ", format(round(maxtemp[i], 1),nsmal=1), "C ",
                                      " Mean: ", format(round(midtemp[i], 1),nsmal=1), "C ",
                                      sep=""), line=1)
      # this plot method is 3 times faster than plot with raster above!!
      # Sys.sleep(0.08) 
      cat('\r',paste("Exporting #",i-1, "of", length(f.start), "frames", sep=" "))
    }, 
    #video.name=paste(f.root,"_animation.mp4",sep=""), # comment out for saveHTML
    #ffmpeg="ffmpeg",                                  # comment out for saveHTML  
    #other.opts = c("-b:v 400k","-pix_fmt yuv420p -crf 18"),   # comment out for saveHTML
    #convert="convert", clean=TRUE,                  # comment out for saveHTML
    img.name = paste(des, "_frame_", sep=""),        # needed for saveHTML
    #imgdir = paste(outputID_Dir,"/img_dir",sep=""), # needed for saveHTML / comment out for saveGIF
    imgdir = "img_dir",
    htmlfile = paste(f.root, ".html",sep=""),    # needed for saveHTML / comment out for saveGIF
    autobrowse = FALSE, navigator=FALSE,         # needed for saveHTML
    description = des)                          # needed for saveHTML
)
cat('\n')
#dev.off()
# Show last plotted image in main plot window
plot.new()
sub.title<-paste(subtitle[i],i)
par(pin=c(4,3),cex.sub=1.2, cex.main=1)
image.plot(d, useRaster=TRUE, bty="n", col=thermal.palette, main=des,
           xlab="", ylab="", xaxt="n", yaxt="n",
           zlim=c(minrangeset,maxrangeset),
           legend.shrink=0.85, legend.cex=0.85, asp=h/w)
title(sub=sub.title, cex.sub=1, adj=0.5, line=3)
title(cex.sub=0.9, adj=0.5, line=1, sub=paste("Min: ", round(mintemp[i], digits=1), "C ",
                                             " Max: ", round(maxtemp[i], digits=1), "C ",
                                              " Mean: ", round(midtemp[i], digits=1), "C ",
                                              sep=""))


#}  # this belongs to the list.files loop



#14. Plot the thermal images - no animation package
#subtitle<-paste(as.character(extract.times)," Frame #: ")
#system.time(
#for(i in 1:length(f.start))
#  {
#    sub.title<-paste(subtitle[i],i)
#    d<-matrix(rev(alltemperature[f.start[i]:f.end[i]]), nrow=w, ncol=h)
#    # populate temperature matrix in reverse in order to flip image before plotting
#    d<-flip.matrix(d)
#          #  ds<-(d-mintemp)/(maxtemp-mintemp) 
#          #par(pin=c(6,4.5)) 
#          #system.time(
#          #  plot(dr,add=FALSE, useRaster=TRUE, axes=FALSE, box=FALSE, xlim=c(0,640), ylim=c(0,480), 
#          #       zlim=c(mintemp,maxtemp), legend.shrink=0.83, legend.width=1.5, col=thermal.palette, 
#          #       main=paste(maintitle[i],i))
#          #          )
#    par(pin=c(4,3),cex.sub=1, cex.main=1)
#    image.plot(d, useRaster=TRUE, col=thermal.palette,  bty="n",
#                 xlab="", ylab="", xaxt="n", yaxt="n", zlim=c(mintemp,maxtemp),
#                 legend.shrink=1, legend.cex=0.85, asp=h/w)
#    title(main=sub.title)
#    # this plot method is 3 times faster than plot with raster above!!
#    Sys.sleep(0.08) 
#  }
#  )