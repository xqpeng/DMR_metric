# DMR_metric
This project provides the service for calculating two metrics (Qn and Ql) for DMR sets predicted by different methods.
In the programs of calculating Qn and Ql, the methylation profiles of samples are needed, which are used to recalculate the methylation levels of CpGs in DMRs. For all the DMRs predicted by different methods, the methylation differences of DMRs are recalculated based on the methylation profiles of samples also. 


There are two programs of preprocessing before calculating Qn and Ql.


Step 1: Merge the methylation report files of samples in each group into one methylation report file     
     Commands: perl merge_group_meth_inf.pl list_methylation_report_files meth_column unmeth_column
      Input: 1. A file containing the methylation report files of samples with group id
          For example: methylation_bedGraph_sample1 A
                methylation_bedGraph_sample2 A
                methylation_bedGraph_sample3 B
                methylation_bedGraph_sample4 B
                ...
          2. the column id to specify the meth counts
          3. the column id to specify the unmeth counts
      Output: merged_methylation_bedGraph_A_B
         
       
Step 2: Format and sort the DMR set predicted by different DMR detection tools
      Commands: perl DMR_format.pl list_DMR_files
      Input: a file containing the DMR files of different methods 
      Output: formated DMR files 
      
Step 3: Calculate Qn, Ql, and Qf under different threshold for each DMR set
      Commands: perl calcluate_QnQl.pl list_formated_DMR_files merged_methylation_bedGraph_A_B
      Input: 1. a file containing the formated DMR files of different methods 
             2. merged_methylation_bedGraph_A_B from step 1.
      Output: 1. a set of DMR files, containing additional information, including the methyation difference and count the CpG number and the length of each DMR
			        2. a file containing Qn, Ql, and Qf of different methods under different threshold 
