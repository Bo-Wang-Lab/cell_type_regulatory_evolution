import pandas as pd
import numpy as np
import gzip
import os
from numba import jit
from multiprocessing import Pool
import argparse

# x is a range array formatted with columns as 
# ID, chrom, start, stop, quantile
# and sorted by quantile
@jit( nopython=True )
def iterOverlapArray( x ):
    # which peaks from this set should we keep
    keep = np.array([False]*x.shape[0])
    # which peaks from this set have already been processed
    seen = np.array([False]*x.shape[0])
    
    for i in range(x.shape[0]):
        if seen[i]:
            continue
        keep[i] = True
        row = x[i,:]
        # Pull out all peaks that overlap this one
        overlap = (x[:,1]==row[1]) &\
                  (x[:,3]>row[2]) &\
                  (x[:,2]<row[3])
        # Ignore the rest (we should be seeing the 
        # 'most important' one first anyway because
        # of the way the array is sorted)
        seen = seen | overlap
    return x[keep,0]

# Read in a list of separate peak files into a single dataframe
# Normalize the score within each file so that we can compare 
#   relative importance of peaks between files
# @param normCol: which column to normalize on - either integer
#   value or one of "score", "signalValue" (case-insensitive)
# @param dropChrs: list of chromosomes to be excluded from analysis
def readBeds( bedlist, gzipped=True, normCol="score", dropChrs=[] ):
    # Check valid score column
    if normCol.lower() == "score":
        normCol = 4
    elif normCol.lower() == "signalvalue":
        normCol = 6
    elif isinstance( normCol, int ) and normCol not in [4,6]:
        print( "Using column number {0}, are you sure you don't want 'score' or 'signalValue'?".format(normCol) )
    else:
        normCol = 4
        print( "normCol must be 'score', 'signalValue', or integer column index, using 'score'" )
    
    # read in the files
    agg = None
    for path in bedlist:
        if gzipped:
            op = gzip.open
        else:
            op = open
        with op( path, "r" ) as f:
            tmp = pd.read_csv( f, sep="\t", header=None )
            tmp = tmp.loc[~tmp.iloc[:,0].isin(dropChrs),:]
            tmp["Quantile"] = tmp.iloc[:,normCol].rank() / len(tmp.index)
            if agg is None:
                agg = tmp
                continue
            agg = pd.concat( (agg,tmp) ).reset_index(drop=True)

    agg = agg.sort_values(by="Quantile",ascending=False).reset_index(drop=True)
    agg.columns = ["Chr","Start","Stop","Name","Score","Strand","Signal","p","q","peak","Quantile"]
    return agg

# Trim a Pandas peak table to have a fixed width for all peaks
# By default, peaks will be centered at their summit
# Peaks which would otherwise over-run the end of a contig
#   are shifted to still be of appropriate length, but be 
#   entirely contained within the sequence
# @param peaks: peak DataFrame formatted as output by 
#   readBeds
# @param targetSize: fixedWidth of peaks
# @param chrSizes: pandas Series whose index is chromosome/contig
#   names and whose values are the corresponding lengths
def resizePeakTable( peaks, targetSize, chrSizes=None ):
    # Try to center everything at the summit
    offset = targetSize // 2
    adj = peaks.copy()
    adj.loc[:,'Start'] = adj.peak - offset
    adj.loc[:,'Stop'] = adj.Start + targetSize
    # Move up peaks that bleed off the start
    tooShort = adj.Start < 0
    adj.loc[tooShort,['Start','Stop']] = adj.loc[tooShort,['Start','Stop']] - adj.loc[tooShort,'Start']
    # Move down peaks that bleed off the end
    if chrSizes is not None:
        print( adj.Stop.values )
        print( chrSizes )
        print( chrSizes.loc[adj.Chr.values].values )
        adj['Over'] = adj.Stop.values - chrSizes.loc[adj.Chr.values].values.flatten()
        tooLong = adj.Over > 0
        adj.loc[tooLong,['Start','Stop']] = adj.loc[tooLong,['Start','Stop']] - adj.loc[tooLong,'Over']
        adj = adj.drop( 'Over', axis=1 )
    return adj

# Perform iterative overlapping of the complete peakset
# Do this as quickly as possible 
# @param peaks: peak DataFrame formatted as output by
#   readBeds
# @param nthreads: how many cores to use - defaults to as
#   many as present on the machine. If 1 or None, use a 
#   single process
def iterativeOverlap( peaks, nthreads=os.cpu_count() ):
    # Convert the index, chromosome, start, stop, and quantile
    # values into a numpy array for faster processing
    asarray = np.zeros((len(peaks.index),5))
    asarray[:,2:] = peaks.iloc[:,[1,2,-1]]
    asarray[:,0] = peaks.index.values

    # assign unique indices to each genome contig
    chrs = peaks.iloc[:,0].unique()
    chrMap = { chrs[x]: x for x in range( peaks.iloc[:,0].nunique() ) } 
    asarray[:,1] = np.array([chrMap[c] for c in peaks.iloc[:,0]])
    
    # run the function once to make sure numba has compiled
    _ = iterOverlapArray( asarray[:1,:] )
    
    # run the full computation, using parallelization if possible
    if nthreads is None or nthreads == 1:
        keep = iterOverlapArray( asarray )
    else:
        # if we're parallelizing, chunk by chromosome since there's never
        # overlap between different contigs
        chunks = [ asarray[asarray[:,1]==c,:] for c in np.unique(asarray[:,1])]
        with Pool( processes=nthreads ) as p:
            keeps = p.map( iterOverlapArray, chunks )
        keep = np.concatenate(keeps)
        
    # Only keep first 10 columns to be a valid narrowPeak file
    # union = peaks.iloc[keep,:-1].sort_values(by=["Chr","Start","Stop"]).reset_index(drop=True)
    union = peaks.iloc[keep,:10].sort_values(by=["Chr","Start","Stop"]).reset_index(drop=True)
    return union

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="Find consensus peaks from one or" + \
                                      " more narrowPeak files using the iterative overlap method" )
    parser.add_argument( "--files", "-i", dest="fnames", type=str, nargs="+", required=True )
    parser.add_argument( "--compare", "-c", dest="normCol", default="score" )
    parser.add_argument( "--threads", "-t", dest="nthreads", default=os.cpu_count(), type=int )
    parser.add_argument( "--exclude", "-e", dest="dropChrs", type=str, nargs="+", default=[] )
    parser.add_argument( "--output", "-o", dest="out", type=str, required=True )
    parser.add_argument( "--is-list", "-L", dest="isList", default=False, action="store_true" )
    parser.add_argument( "--resize", "-r", dest="resize", default=0, type=int )
    parser.add_argument( "--chrom-sizes", "-z", dest="sizes", default="", type=str ) 
    
    args = parser.parse_args()
    
    normCol = args.normCol
    nthreads = args.nthreads
    dropChrs = args.dropChrs
    out = args.out
    isList = args.isList
    
    files = args.fnames
    if isList:
        listName = files[0]
        with open(listName, "r") as f:
            files = [str(l).strip() for l in f.readlines()]
    files = list(filter( lambda x: len(x) > 0, files ))
    gzipped = files[0].endswith(".gz")
    
    print( "Importing peaks from {0} files".format(len(files)) )
    
    agg = readBeds( files, dropChrs=dropChrs, gzipped=gzipped, normCol=normCol )
    
    print( "{0} total (potentially-overlapping) peaks identified".format(len(agg.index)) )
    
    resize = args.resize
    sizePath = args.sizes
    sizes = None
    if len(sizePath) > 0:
        sizes = pd.read_table( sizePath, header=None, index_col=0 )
    if resize > 0:
        print( "Resizing peaks to {0} bases, centered on the summit".format(resize) )
        agg = resizePeakTable( agg, resize, sizes )

    print( "Collapsing overlaping peaks" )
    union = iterativeOverlap( agg, nthreads=nthreads )
    
    print( "{0} non-overlapping peaks remaining after iterative overlap removal".format(len(union.index)) )
    
    union.to_csv( out, sep="\t", index=False, header=False )
