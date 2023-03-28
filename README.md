# splicingPatternMotifAnalysis
The script for using meme to explore motif enrichment in different splicing patterns. Needed to be used in tandem with shiba

# Usage:
The script will try to access txt files under results which is produced by shiba. Typically there are 5 files which correspond to 5
splicing patterns that we are interested in. Then this program will call meme to run the analised result. The final output will be
in memeresult folder. 

# Caution:
The program will return the shell (between some time and) immediately upon launch, the process is still running! But the shell will be
returned. This is the result of the program trying to multithreading everything, which detaches the stdio from the shell. While it
might have returned, the program is still running in the background, please do not close the shell otherwise the detached program
will be killed. To get some understanding of the computation, use htop or top and the users will be able to locate the active computing
process until all of them are finished. Since all meme instances are launched for all processed sampled, the host might run out of memory.
Tested on 90GB server, peak memory usage 30 GB. A node with 50GB memory should suffice.
