How to run TrackIsolation code
==============================

Running method of all codes is the same as submtting signal window code. But in PiXTRK code, BDT egamma ID of HGCAL is added. So someone who wants to run codes, he/she must use these variables: <br>


For signal sample: <br>
<blockquote>
<p> InputDir = "/xrootd/store/user/jhkim/L1Pixel/SE_PU200_HGCAL_EGID" </p>
<p> number_of_cores = 269 </p>
<p> fr = open('./inputlist_SE_PU200_HGCAL_EGID.txt', 'r') </p>
<blockquote>


For minbias sample: <br>
<blockquote>
<p> InputDir = "/xrootd/store/user/jhkim/L1Pixel/MinBias_PU200_HGCAL_EGID" </p>
<p> number_of_cores = 621 </p>
<p> fr = open('./inputlist_MinBias_PU200_HGCAL_EGID.txt', 'r') </p>
<blockquote>
