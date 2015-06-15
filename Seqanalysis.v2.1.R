########### TABLE OF CONTENTS ############################
#0. Define all libraries and functions
#1. File Handling - manually enter file name as variable f
#2. Extract meta-tags from thermal vid file
#3. File information for binary read in
#4. Read binary Data in (consider using n=c())???
#5. Find everywhere in vector where resolution info is: i.e. 640x480, 320x240
#6. Fill the array 'header' with values.  Create an index variable for the frame data
#7. Extract the time values from the headers
#8. Convert raw data held in alldata to temperature
#9. Image variability analysis to estimate activity difference frame by frame.
#10. Convert activity data into cumulative sums 
#11. Calculate slope of the cumulative function every N samples
#12. Put final activity data into csv files
#13. Plot the thermal images.
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# Version History
# Version 0-1 disorganised assortment, trials associated with opening up raw files 
# Version 1-2  Figured out the math behind raw to temperature conversion.  
# Version 2. Consolidated and annotated, added table of contents.
#     Added video output. Added frame-frame activity analysis and subsampling of activity
#     Added enhancement to frame-indexing, removing for-loop slowness by implementing unlist(map())

#0. Define functions & required libraries ####
library("lattice")
library("fields")
library("animation")
library("Thermimage")
library("raster")

# !!! Read this section before starting  !!!!
# Requires: in addition to the libraries listed above:
# Need to have exiftool installed in your OS's system folder or equivalent
# For Windows: Copy the exiftool.exe file to the C:/windows folder and bestow administrator
# priviledges
# For OSX: Simply install the .dmg file.  It will be installed automatically into the 'system'
# For UNIX/Server: download the PERL script and place in home folder: "~/src/exiftool"
# http://www.sno.phy.queensu.ca/~phil/exiftool/

#--- Set the colour palette to be used  --- --- --- --- --- --- --- --- --- --- --- --- 
# for use with image.plot and other rasterised thermal images
thermal.palette<-palette.choose("ironbow")  # can choose form "flir","ironbow"...need to add others
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#1. File Handling - manually enter file name as variable f ####
# define the file output names and locations
f<-NULL
mainDir<-"~/Desktop/thermalvids/"
outputDir<-"/Documents/R/Scripts and Analysis/Analysis/Galapagos/activity/"
getwd() #show working directory
setwd(mainDir)              

txt.outDir<-paste(outputDir,"txt_output",sep="")
# data arising from this script will be saved here
dir.create(txt.outDir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
vid.outDir<-paste(outputDir,"vid_output/",sep="")
# video/image outputs from this script with be saved here
dir.create(vid.outDir, showWarnings = TRUE, recursive = FALSE, mode = "0777")

# l.files<-list.files(mainDir,all.files=FALSE,full.names=FALSE,include.dirs=TRUE)

#f<-"Sandpiper0001.SEQ"
f<-"MAG13-Night0001.SEQ"
#f<-"Metronome40bpm0001.SEQ"

f.root<-substr(f,1,nchar(f)-4)             # text output file name root, without .rtv/seq
f.raw.activity<-paste(data.outDir, f.root,"_raw.activity.csv", sep="")
f.resamp.activity<-paste(data.outDir, f.root,"_resamp.activity.csv", sep="")
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#2. Extract meta-tags from thermal vid file ####
# source: http://timelyportfolio.github.io/rCharts_catcorrjs/exif/
# see also here for converting thermal image values
# http://u88.n24.queensu.ca/exiftool/forum/index.php?topic=4898.45
# accessing exiftool from system command line
# Decipher Camera Meta Data information 
# Need to have exiftool installed in your OS's system folder or equivalent
# http://www.sno.phy.queensu.ca/~phil/exiftool/
camvals<-"-flir -*Emissivity -*Original -*Date -*Planck* -*Distance -*Temperature* -*Transmission -*Humidity -*Height -*Width -*Model* -*Median -*Range"
if (Sys.info()["sysname"]=="Darwin")
  {
  info<-system(paste("exiftool", f, camvals),intern=TRUE)
  } 
# if trouble calling exiftool, try using "DIR/exiftool":
if (Sys.info()["sysname"]=="Linux")
 {
  info<-system(paste("~/src/exiftool/exiftool", f, camvals),intern=TRUE)
 }
img.df<-read.fwf(textConnection(info), widths=c(32,1,1,60), stringsAsFactors=FALSE, header=FALSE, fill=FALSE)
img.df<-img.df[,-c(2,3)]
values<-gsub("[^0-9: .-]","", img.df[,2])
suffixes<-gsub("([0-9.])","",img.df[,2])
img.df<-img.df[,-2]
img.df<-data.frame(img.df,values,suffixes)
img.df<-as.matrix(img.df)

# Set variables for calculation of temperature values from raw A/D sensor data 
op<-options(digits.secs=3)
w<-NULL;h<-NULL
Emissivity<-  as.numeric(img.df[1,2])       # Object Emissivity - should be ~0.95 or 0.96
dateOriginal<-as.character(img.df[2,2])
dateModif<-   as.character(img.df[3,2])
GMT<-         substr(img.df[2,2],nchar(img.df[2,2])-5,nchar(img.df[2,2]))
PlanckR1<-    as.numeric(img.df[6,2])       # Planck R1 constant for camera  21106.77
PlanckB<-     as.numeric(img.df[7,2])       # Planck B constant for camera  1501
PlanckF<-     as.numeric(img.df[8,2])       # Planck F constant for camera  1
PlanckO<-     as.numeric(img.df[9,2])       # Planck O constant for camera  -7340
PlanckR2<-    as.numeric(img.df[10,2])      # Planck R2 constant for camera 0.012545258
OD<-          as.numeric(img.df[11,2])      # object distance in metres
FD<-          as.numeric(img.df[12,2])      # focus distance in metres
ReflT<-       as.numeric(img.df[13,2])      # Reflected apparent temperature
AtmosT<-      as.numeric(img.df[14,2])      # Atmospheric temperature (for loss across distance)
IRWinT<-      as.numeric(img.df[15,2])      # IR Window Temperature  20
IRWinTran<-   as.numeric(img.df[18,2])      # IR Window transparency 1
RH<-          as.numeric(img.df[19,2])      # Relative Humidity  50
h<-           as.numeric(img.df[20,2])      # sensor height (i.e. image height)
w<-           as.numeric(img.df[21,2])      # sensory width (i.e. image width)
rawmedian<-   as.numeric(img.df[25,2])     # raw median value (first frame?)
rawrange<-    as.numeric(img.df[26,2])     # raw range value (first frame?)
rawmax<-      rawmedian+rawrange/2
rawmin<-      rawmedian-rawrange/2
                         
op<-options(digits.secs=3)
dateOriginal<-substr(dateOriginal,1,nchar(dateOriginal)-6)
datestring<-gsub(":","-",substr(dateOriginal,1,10))
timestring<-substr(dateOriginal,12,23)
dateOriginal<-paste(datestring,timestring)
dateOriginal<-strptime(dateOriginal,"%Y-%m-%d %H:%M:%OS")

# Set thermal image frame parameters here w = width, h = height
# you must write these values in or the script will not work!
if(is.null(w)) w<-640 
if(is.null(h)) h<-480
l<-w*h
footerlen<-52
# not all seq files have the same footer duration, although I have not checked for certain
# this might be a camera specific feature.
#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#3. File information for binary read in ####
finfo <- file.info(f)
# finfo$size  # how many bytes in the file
byte.length<-2  # how many bytes make up an element
no.elements<-finfo$size/byte.length 
# this should tell you how many total elements make up the inputted data
op<-options(digits.secs=3)
save.time<-strptime(finfo$mtime,"%Y-%m-%d %H:%M:%OS")
# time at which the file was saved.  this almost corresponds to the last frame of the file
if(dateOriginal<save.time) save.time<-dateOriginal

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#4. Read initial 100,000,000 bytes binary Data in  ####
to.read <- file(f, "rb") # set to.read file.  rb means read binary
system.time(alldata<-readBin(to.read, integer(), n=100000000, 
                     size=byte.length, endian = "little", signed=FALSE))
no.frames<-round(length(alldata)/l)
# approx number of frames in video file, can't be certain since some headers might be different sizes
# this value might be an over-estimate
close(to.read)
# alldata[1:1000]
# on first read, print out the first 1000 elements of alldata
# you can tell the actual thermal data, since every data point for a broad swath is a 4 to 5 digit number
# that could range up to 15 or 16 bits (2^15) for the raw data the thermal camera records
# this is such that a value of 0 would be the lowest measurement and a value of 2^16 the highest value
# the thermal camera can detect 

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#5. Find everywhere in vector where resolution info is: i.e. 640x480, 320x240 ####
# note: user must set l and w at the top of this script
fid<-c(w,h)             
# this is the look-up sequence, which will repeat throughout the file in a predicable fashion
# trouble with very large files, so try a small portion and predict the remaining 
# wh.locates

if(length(alldata)>=10000000)
  {
  system.time(wh.locate<-locate.fid(fid,alldata,long=TRUE))
  # try wh.locate on a small chunk of all data
  diff.wh.locate<-diff(wh.locate) 
  # difference calc should yield a repeating pattern of header then frame
  gaps<-unique(diff.wh.locate)
  # if the pattern is simple, this value should be 2
  no.unique.locates<-length(gaps) 
  # reconstruct wh.locate from scrap, starting only with wh.locate[1]
  repeats<-trunc(finfo$size/2/(l+wh.locate[1]+gaps[1]))
  # how many repeats required to create fill whole file
  wh.locate<-cumsum(as.numeric((c(wh.locate[1],rep(gaps,repeats)))))
  # cumulative sum up the 1st locate and repeate gaps
  wh.locate<-(wh.locate[-c(which(wh.locate>(finfo$size/byte.length)))])
  # remove any locates that go beyond the length of the data file after import
  
}

if(length(alldata)<10000000)
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

#header.l<-diff.wh.locate[length(diff.wh.locate)]
header.l<-as.integer(rev(wh.locate)[1]-rev(wh.locate)[2])
# the approximate length of the header, as defined by the last element in the diff.wh.locate vector
# all thermal vids appear to end with a 'header' with two resolution statements spaced by this amount
# the frame data stream commences 15 elements after the resolution info
res2fram<-15
# Below define the start indices for the headers (h.start) and frames (f.start)
# check if the first location of the resolution info is a small value or not
# .SEQ files have two instances where resolution is recorded, .FCF files only have it once toward
# the end of the header

if(wh.locate[1]<header.l & no.unique.locates==2)  # .SEQ files appear to be formatted this way  
{
  h.start<-wh.locate-header.l
  h.start<-h.start[seq(2,length(h.start),2)]
  h.end<-h.start+header.l+res2fram-1
  f.start<-wh.locate+res2fram
  f.start<-f.start[seq(2,length(f.start),2)]
  f.end<-f.start+l-1
} else if(wh.locate[1]>header.l & no.unique.locates==2) 
  # some .fcf formatted this way
{
  h.start<-wh.locate-header.l
  h.start<-h.start[seq(2,length(h.start),2)]
  h.end<-h.start+header.l+res2fram-1
  f.start<-wh.locate+res2fram
  f.start<-f.start[seq(2,length(f.start),2)]
  f.end<-f.start+l-1
} else if (wh.locate[1]>=header.l & no.unique.locates>2) 
  # other .fcf files formatted this way
{
  wh.locate<-wh.locate[-2]
  h.start<-wh.locate-header.l
  h.start<-h.start[seq(1,length(h.start),2)]
  h.end<-h.start+header.l+res2fram-1
  f.start<-wh.locate+res2fram
  f.start<-f.start[seq(1,length(f.start),2)]
  f.end<-f.start+l-1
} 

# each f.start should correspond to the start of a frame in the video file
# from my impression, the location of the header is 15 elements in front of the first pixel value
# above, res2fram is set to be 15

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#6a. Read header information in with readBin and seek ####
# since we have defined the h.starts and h.ends above, in principle this can be used
# to open the "to.read" variable at only those locations
temp<-NULL
timeplace<-NULL
to.read <- file(f, "rb") # set to.read file.  rb means read binary
for (i in 1:length(h.start))
{
  seek(to.read,where=(h.start[i]+448)*byte.length,origin="start")
  temp<-readBin(to.read, integer(), n=3, 
              size=byte.length, endian = "little", signed=FALSE)
  timeplace<-c(timeplace,temp)
  cat('\r', paste("Timeplace # ", i, sep=""))
}
close(to.read)
timeplace<-unname(matrix(timeplace, nrow=3, byrow=FALSE))

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#6b. Read/seek and fill the frame with data.   ####
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
  cat('\r', paste("Frame # ", i, sep=""))
}
)
close(to.read)
cat('\n')
cat(paste("The following frames were flagged for errors:"))
flagged.fram
flagged.fram.index<-rep(0,length(f.start))
flagged.fram.index[flagged.fram]<-1

#header<-NULL
# generate header data
#  header<-alldata[h.start[1]:h.end[1]]
#  p.ind<-seq(1,h.start[1]-1,1)
#  # define the initial part of all data prior to the normal header start
#  #h.ind<-h.start[1]:h.end[1]
# system.time(for (i in 2:length(h.start))
#  {
#    j<-i-1 
#    h.temp<-alldata[h.start[i]:h.end[i]]
#    #h.ind.temp<-h.start[i]:h.end[i]
#    p.ind.temp<-seq(f.end[j]+1,h.start[i]-1,1)
#    header<-cbind(header,h.temp)
#    #h.ind<-c(h.ind,h.ind.temp)
#    p.ind<-c(p.ind,p.ind.temp)
#    
#  }
# )
#  # check to make sure that the length of the indices divded by the length of the header 
#  # equals an even integer equal to the number of frames in the file

# a fast way of looping through both vectors (f.start and f.end) and creating a sequence with
# the ':' command
# in principle, the loop above is now redundant and only h.ind and f.ind below is important
#system.time(h.ind<-unlist(Map(':',h.start, h.end)))
#system.time(header<-matrix(alldata[h.ind],nrow=h.end[1]-h.start[1]+1))

#system.time(f.ind<-as.numeric(unlist(Map(':',f.start, f.end))))
# consider a faster way?  Maybe:

#system.time(f.ind<-as.numeric(unlist(Map(':',f.start,f.start+l-1))))
# have removed f.ind - killng memory


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#7. Extract the time values from the headers ####
# Technically, i could re-load the data 1 byte at a time and find the date stamps and convert as I do above
# where I found dateOriginal, but this is an alternative
# Once you have the headers defined, do a subtraction of the first two headers to find
# the location of the millisecond marker place

# framediffcheck<-header[,2]-header[,1]
# framediff<-apply(header,1,diff)
# if all goes well, only one element will show a difference
# for my test .SEQ files, this turns out to be the 452nd element
# or more accurately, 137 elements BEFORE the frame starts
# thus: 452 contains msec, 450 contains sec, 451 contains day info
# so far, it is the same for my test .fcf file as well
# these 3 values have subsequently been extracted and placed into a new variable
# called timeplace, where
# timeplace[2] = dayplace, timeplace[1]=secplace, timeplace[3]=msecplace
#msecplace<-452
#secplace<-msecplace-2
#dayplace<-msecplace-1

time.msec<-as.vector(timeplace[3,])
time.sec<-as.vector(timeplace[1,])
#if(time.sec[2]-time.sec[1]>1) time.sec[1]<-time.sec[2]
# may or may not require this if statement any longer, but sometimes the first header
# lacks information on the # of elapsed minutes
time.sec<-time.sec
f.times<-time.sec+time.msec/1000
f.times<-f.times-f.times[1]
# f.times should provide a vector corresponding to the time in elapsed seconds for each frame
# if all worked out properly, that is!
options(digits.secs=3)
est.times<-save.time+f.times
#f.times<-f.times-f.times[1]
# real times in Date:Time format based on finfo$mtime which should correspond to the time
# the file was saved.  Only problem with this is that if you edit or re-save a .seq file, then
# this will not work

# I think I found the way to convert these time values a different way
# help from http://www.silisoftware.com/tools/date.php
# http://www.sandersonforensics.com/forum/content.php?131-A-brief-history-of-time-stamps

hex.time.day<-as.hexmode(timeplace[2,])
hex.time.sec<-as.hexmode((timeplace[1,]))
time.char<-paste("0x",hex.time.day,hex.time.sec,sep="")
time.num<-as.numeric(time.char)+timeplace[3,]/1000
extract.times<-as.POSIXct(time.num, origin="1970-01-01")
rownames(extract.times)<-NULL

#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#8. Convert raw data held in alldata variable to temperature ####
# use parameters relevant to experiment.  Distance 0.5 m, Feather/body emissivity 0.96
# optimal window transmittance 0.96, temperatures set to 20C, RH set to 50%
# Technically the RTemp, ATemp and IRWTemp should be the incubator temperature
templookup<-raw2temp(1:65535,0.96,0.4,20,20,20,0.96,50,PlanckR1,
                     PlanckB,PlanckF,PlanckO,PlanckR2)
# create a lookup vector, might be faster than calculating?

#alltemperature<-raw2temp(alldata,0.96,0.4,20,20,20,0.96,50,PlanckR1,PlanckB,PlanckF,PlanckO,PlanckR2)
# lookup temperature from the lookup vector...eliminate calculating
rangetemp<-NULL
rangetemp[1:length(f.start)*2]<-NA
j<-seq(1,length(f.start)*2,2)
system.time(
  for (i in 1:length(f.start))  
{
  rng.temp<-range(templookup[alldata[1:(l-1),i]])
  rangetemp[j[i]:(j[i]+1)]<-rng.temp
  cat('\r',paste("Calculating Frame #", i, "min:", round(rng.temp[1],digits=1), "max:", 
            round(rng.temp[2], digits=1), "and mid:", round(sum(rng.temp)/2,digits=1), "     ", sep=" "))
}
) 
mintemp<-rangetemp[c(seq(1,length(rangetemp),2))]
maxtemp<-rangetemp[c(seq(2,length(rangetemp),2))]
midtemp<-(maxtemp+mintemp)/2
allmintemp<-median(mintemp, na.rm=TRUE)
allmaxtemp<-median(maxtemp, na.rm=TRUE)
allmidtemp<-median(midtemp, na.rm=TRUE)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

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
  cat('\r',paste("Calculating differences between Frame# ",i-1, " and ", i,
                 which(d==NA), sep=""))
  #Sys.sleep(0.01)
  #if(i>1 & sd(d)==0)  
    # break out of for loop the first instance of a frame freeze occurs (common in Mikron
    # when there is a memory problem during data recording and all images are identical)
  #{
  #  break
  #}
  av<-mean(d)           # average of the |inter-frame difference|
  sde<-sd(d)            # sd of the |inter-frame difference|
  avg<-c(avg,av)        # bind the current av to previous avg to create a vector of avg abs(diff) values
  sdev<-c(sdev,sde)     # bind the current sde to previous sdev to create a vector of sd of abs(diff) values
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#10. Convert activity data into cumulative sums ####
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
matplot(df.time, cumsum.act, type = "l",  lty=c(1,1), lwd=1, xlab="Time (s)", ylab="Cumulative Sum",
        ylim=c(min(cumsum.act), max(cumsum.act)), xlim=c(min(df.time),max(df.time)))
legend("topleft", bty = "n", c("Avg", "SD"), lty=c(1,1), col=c("black", "red"))
# plot the cumulative function 

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#11. Calculate slope of the cumulative function every N samples ####
samples<-60 # set this value to the number of sample points over which to calculate avg(slope) data 
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
samp.sec<-meanEveryN(dfsec,samples)
#samp.sec<-(samp.sec[,-1])
samp.min<-samp.sec/60
matplot(samp.min, slp[,2:3], type = "l",  lty=c(1,1), lwd=1, xlab="Time (min)", ylab="Slope of cumulative sum",
        ylim=c(min(slp[,2:3]), max(slp[,2:3])), xlim=c(min(dfmin),max(dfmin)))
legend("topleft", bty = "n", c("Avg", "SD"), lty=c(1,1), col=c("black", "red"))

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#12. Put final activity data into csv files ####
nf<-length(f.start) # quick count of the number of frames to be exported 
# although this value (nf) is the total frame #.  for the analysis of frame difference analysis: nf-1
activity.output<-data.frame(2:nf, flagged.fram.index[2:nf], extract.times[2:nf], est.times[2:nf], f.times[2:nf], activity, mintemp[2:nf],
                            maxtemp[2:nf], midtemp[2:nf])
colnames(activity.output)<-c("Frame #", "Flagged Frame", "Ext. Times", "Est. Times", "Elapsed Time(s)", "Average", "SD", "Min Temp", "Max Temp", "Mid Temp")
resam.output<-data.frame(samp.times,samp.min,samp.sec,slp)
colnames(resam.output)<-c("Samp Times", "Elapsed Time (min)", "Elapsed Time(s)","Sample #", "Slope Avg","Slope SD")
write.csv(activity.output, f.raw.activity)
write.csv(resam.output, f.resamp.activity)
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#13. Plot the thermal images. ####
# This requires the "animation" package, which is not available in ready form on CRAN. 
# Therefore, you must download the .tar file to your computer and install from there.
# install.packages("~/Downloads/animation_2.3.tar", repos= NULL, type="mac.binary.mavericks")  
# library("animation") need to download seewave tar from website:
# http://www.renevolution.com/how-to-install-ffmpeg-on-mac-os-x/
# Save image sequence as video, with image number of and time. saveVideo
#Or use saveHTML function
#Then Open images in R and take thermometrics. (not implemented yet)
des = c(f.root)
setwd(vid.outDir)
subtitle<-paste(as.character(extract.times)," Frame #: ")
system.time(
  #saveHTML(  
  #saveGIF(  
  saveVideo( 
  for(i in 1:length(f.start))
    {
      #ani.options(interval = 0.08)
      ani.options(interval = 1/Fs, ani.width=720)
      sub.title<-paste(subtitle[i],i)
      #fend<-f.start[i]+l-1
      #t<-raw2temp(alldata[,i],0.96,0.4,20,20,20,0.96,50,PlanckR1,
       #                       PlanckB,PlanckF,PlanckO,PlanckR2)
      #d<-matrix(rev(t), nrow=w, ncol=h)
      d<-matrix(rev(templookup[alldata[,i]]), nrow=w, ncol=h)
        # use the templookup variable, to reduce the requirement for alltemperature
        # populate temperature matrix in reverse in order to flip image before plotting 
      d<-flip.matrix(d)
      #dr<-raster(rotate90.matrix(d), xmn=0, xmx=640, ymn=0, ymx=480) # convert matrix to raster
      #  ds<-(d-mintemp)/(maxtemp-mintemp) 
      #plot(dr,add=FALSE, useRaster=TRUE, box=FALSE, axes=FALSE, main=paste(maintitle[i],i), 
      #     xlim=c(0,640), ylim=c(0,480), zlim=c(mintemp,maxtemp), 
      #     legend.shrink=0.8, legend.width=1, col=thermal.palette)
      par(pin=c(6,4.5),cex.sub=1.5, cex.main=1.5)
      image.plot(d, useRaster=TRUE, bty="n", col=thermal.palette,
                 xlab="", ylab="", xaxt="n", yaxt="n", 
                 zlim=c(allmintemp,allmaxtemp),
                 legend.shrink=0.85, legend.cex=0.85, asp=h/w)
      title(main=des, sub=sub.title)
      title(cex=1, sub=paste("Min: ", round(mintemp[i], digits=1), "C ",
                             " Max: ", round(maxtemp[i], digits=1), "C ",
                             " Mean: ", round(midtemp[i], digits=1), "C ",
                             sep=""), line=1)
      # this plot method is 3 times faster than plot with raster above!!
      # Sys.sleep(0.08) 
      cat('\r',paste("Exporting #",i-1, "of", length(f.start), "frames", sep=" "))
    }, 
    video.name=paste(f.root,"_animation.mp4",sep=""), # comment out for saveHTML
    ffmpeg="ffmpeg",                                  # comment out for saveHTML  
    other.opts = c("-b:v 400k","-pix_fmt yuv420p -crf 18"),   # comment out for saveHTML
    #convert="convert", clean=TRUE,                   # comment out for saveHTML
    img.name = "Frame_",                        # needed for saveHTML
    #imgdir = paste(vid.outDir, f.root, ".img_dir",sep=""), # needed for saveHTML / comment out for saveGIF
    #htmlfile = paste(f.root, ".html",sep=""),  # needed for saveHTML / comment out for saveGIF
    #autobrowse = TRUE, navigator=FALSE,        # needed for saveHTML
    description = des)                          # needed for saveHTML
)
dev.off()
graphics.off()
# Show last plotted image in main plot window
par(pin=c(6,4.5),cex.sub=1.5, cex.main=1.5)
image.plot(d, useRaster=TRUE, bty="n", col=thermal.palette, main=des,
           sub=sub.title, line=1, xlab="", ylab="", xaxt="n", yaxt="n", 
           zlim=c(allmintemp,allmaxtemp),
           legend.shrink=0.85, legend.cex=0.85, asp=h/w)
title(cex=1, sub=paste("Min: ", round(mintemp[i], digits=1), "C ",
                " Max: ", round(maxtemp[i], digits=1), "C ",
                " Mean: ", round(midtemp[i], digits=1), "C ",
                sep=""), line=3)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


#14. Plot the thermal images - non gif animation. ####
# This requires the "animation" package, which is not available in ready form on CRAN. 
# Therefore, you must download the .tar file to your computer and install from there.
# install.packages("~/Downloads/animation_2.3.tar", repos= NULL, type="mac.binary.mavericks")  
# library("animation") need to download seewave tar from website, but also 
#Save image sequence as video, with image number of and time. saveVideo
#Or use saveHTML function
#Then Open images in R and take thermometrics. (not implemented yet)
#des = c(f.root)
#setwd(vid.outDir)
#maintitle<-paste(as.character(extract.times)," Frame #: ")
#saveHTML(  
#  for(i in 1:length(f.start))
#  {
#    main.title<-paste(maintitle[i],i)
#    ani.options(interval = 0.08)
#    fend<-f.start[i]+l-1
#    d<-matrix(rev(templookup[alldata[f.start[i]:f.end[i]]]), nrow=w, ncol=h)
#       # use the templookup variable, to reduce the requirement for alltemperature
#       # populate temperature matrix in reverse in order to flip image before plotting
#    d<-matrix(rev(alltemperature[f.start[i]:f.end[i]]), nrow=w, ncol=h)  d<-matrix(rev(alltemperature[f.start[i]:fend), nrow=w, ncol=h)
#    d<-rotate90.matrix(d)
#    d<-flip.matrix(d)
#    # raster function comes from library("raster")
#    dr<-raster(d, xmn=0, xmx=640, ymn=0, ymx=480) # convert matrix to raster
#    #dev.off()
#    par(pin=c(6,4.5))   # set size of image 5 inch x 3.75 inch
#    plot(dr,add=FALSE, useRaster=TRUE, axes=FALSE, box=FALSE, xlim=c(0,640), ylim=c(0,480), zlim=c(mintemp,maxtemp),
#         legend.shrink=0.83, legend.width=1.5, col=thermal.palette, main=main.title)
#    ## two comments hashes refer to statements that worked together with a previous iteration
#    ##image1<-as.matrix(ds)
#    #plot(c(0, w), c(0, h), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", main=maintitle, col=rainbow(99, 0, 1)) #set plotting area
#    #rasterImage(image1, 0, 0, w, h, interpolate = TRUE) #plot image
#    ##image(flip(image1), useRaster=TRUE, col=thermal.palette, main=maintitle[i], xlab="", ylab="", xaxt="n", yaxt="n")
#  }, 
#  img.name = "frame_", imgdir = "img_dir", htmlfile = paste(f.root, ".html",sep=""), 
#  autobrowse = TRUE, 
#  description = des)
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


