How to submit Signal window code
================================

Before submitting, you have to open batch\_script.py and modify code to fit your own account or machine. <br>
Especially, **line 9**


1. If you want to use samples at this path: /xrootd/store/user/jhkim/L1Pixel/SE\_PU200 <br>
   You have to match these variables: <br>
   <blockquote>
    <p> InputDir = "/xrootd/store/user/jhkim/L1Pixel/SE_PU200"</p>
    <p> number_of_cores = 225</p>
    <p> fr = open('./inputlist_SE_PU200.txt', 'r')</p>
   </blockquote> 

2. If you want to use samples at this path: /xrootd/store/user/jhkim/L1Pixel/SE\_PU200\_HGCAL\_EGID <br>
   You have to match these variables: <br>
   <blockquote>
    <p> InputDir = "/xrootd/store/user/jhkim/L1Pixel/SE_PU200_HGCAL_EGID"</p>
    <p> number_of_cores = 269</p>
    <p> fr = open('./inputlist_SE_PU200_HGCAL_EGID.txt', 'r')</p>
   </blockquote>
