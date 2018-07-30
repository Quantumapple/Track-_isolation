How to submit Signal window code
================================

Before submitting, you have to open batch\_script.py and modify code to fit your own account or machine. <br>
Especially, **line 9**  <br>
   <blockquote>
    <p> config+='cd /cms/ldap_home/jongho/CMSSW/CMSSW_9_3_7/src \n'</p>
   </blockquote> 


Then, open submit\_control.py to edit the code.

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

3. Depend on your account, change 'workspace' path to your own directory.
   <blockquote>
    <p> workspace = "/cms/ldap_home/jongho/L1PixelTrigger/SignalWindow/Results/"</p>
   </blockquote> 

## Run command
After changing all variables, then you put:
   <blockquote>
    <p> python submit_control.py</p>
   </blockquote>
