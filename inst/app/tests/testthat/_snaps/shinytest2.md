# {shinytest2} recording: Single-Cell Tab

    Code
      init_files
    Output
       [1] "test_data_dir/.gs_dir"                         
       [2] "test_data_dir/.indices_dir"                    
       [3] "test_data_dir/.pert_query_dir"                 
       [4] "test_data_dir/.pert_signature_dir"             
       [5] "test_data_dir/.tx2gene_dir"                    
       [6] "test_data_dir/test_user"                       
       [7] "test_data_dir/test_user/default"               
       [8] "test_data_dir/test_user/default/bulk"          
       [9] "test_data_dir/test_user/default/custom_queries"
      [10] "test_data_dir/test_user/default/single-cell"   
      [11] "test_data_dir/test_user/prev_project.qs"       

---

    Code
      dataset_files
    Output
       [1] "test_data_dir/.tx2gene_dir/hsapiens_tx2gene.qs"                           
       [2] "test_data_dir/test_user/default/single-cell/mock_10x"                     
       [3] "test_data_dir/test_user/default/single-cell/mock_10x/counts.qs"           
       [4] "test_data_dir/test_user/default/single-cell/mock_10x/dgclogs.qs"          
       [5] "test_data_dir/test_user/default/single-cell/mock_10x/founder.qs"          
       [6] "test_data_dir/test_user/default/single-cell/mock_10x/resoln.qs"           
       [7] "test_data_dir/test_user/default/single-cell/mock_10x/shell.qs"            
       [8] "test_data_dir/test_user/default/single-cell/mock_10x/snn1"                
       [9] "test_data_dir/test_user/default/single-cell/mock_10x/snn1/annot.qs"       
      [10] "test_data_dir/test_user/default/single-cell/mock_10x/snn1/applied.qs"     
      [11] "test_data_dir/test_user/default/single-cell/mock_10x/snn1/clusters.qs"    
      [12] "test_data_dir/test_user/default/single-cell/mock_10x/snn1/scseq_sample.qs"
      [13] "test_data_dir/test_user/default/single-cell/mock_10x/snn1/summed.qs"      
      [14] "test_data_dir/test_user/default/single-cell/mock_10x/snn_graph.qs"        
      [15] "test_data_dir/test_user/default/single-cell/mock_10x/species.qs"          
      [16] "test_data_dir/test_user/default/single-cell/mock_10x/tlogs.tenx"          

---

    Code
      contrast_choices
    Output
        test ctrl      name  value testColor ctrlColor            title
      1    1  all CD14 Mono      1   #CDCC5D     white CD14 Mono vs all
      2    1    2 CD14 Mono 1-vs-2   #CDCC5D   #ED665D   CD14 Mono vs 2
      3    1    3 CD14 Mono 1-vs-3   #CDCC5D   #ED97CA   CD14 Mono vs 3
      4    1    4 CD14 Mono 1-vs-4   #CDCC5D   #729ECE   CD14 Mono vs 4
      5    1    5 CD14 Mono 1-vs-5   #CDCC5D   #FF9E4A   CD14 Mono vs 5

---

    Code
      saved_metric_files
    Output
      [1] "test_data_dir/test_user/default/single-cell/mock_10x/snn1/saved_metrics.qs"

---

    Code
      change_resoln_files
    Output
      [1] "test_data_dir/test_user/default/single-cell/mock_10x/snn2"                
      [2] "test_data_dir/test_user/default/single-cell/mock_10x/snn2/annot.qs"       
      [3] "test_data_dir/test_user/default/single-cell/mock_10x/snn2/applied.qs"     
      [4] "test_data_dir/test_user/default/single-cell/mock_10x/snn2/clusters.qs"    
      [5] "test_data_dir/test_user/default/single-cell/mock_10x/snn2/scseq_sample.qs"
      [6] "test_data_dir/test_user/default/single-cell/mock_10x/snn2/summed.qs"      

---

    Code
      integrated_dataset_name
    Output
      [1] "mock_10x_integrated_harmony"

---

    Code
      integrated_files
    Output
       [1] "test_data_dir/test_user/default/single-cell/mock_10x_integrated_harmony"                     
       [2] "test_data_dir/test_user/default/single-cell/mock_10x_integrated_harmony/args.json"           
       [3] "test_data_dir/test_user/default/single-cell/mock_10x_integrated_harmony/counts.qs"           
       [4] "test_data_dir/test_user/default/single-cell/mock_10x_integrated_harmony/dgclogs.qs"          
       [5] "test_data_dir/test_user/default/single-cell/mock_10x_integrated_harmony/founder.qs"          
       [6] "test_data_dir/test_user/default/single-cell/mock_10x_integrated_harmony/pairs.qs"            
       [7] "test_data_dir/test_user/default/single-cell/mock_10x_integrated_harmony/resoln.qs"           
       [8] "test_data_dir/test_user/default/single-cell/mock_10x_integrated_harmony/shell.qs"            
       [9] "test_data_dir/test_user/default/single-cell/mock_10x_integrated_harmony/snn1"                
      [10] "test_data_dir/test_user/default/single-cell/mock_10x_integrated_harmony/snn1/annot.qs"       
      [11] "test_data_dir/test_user/default/single-cell/mock_10x_integrated_harmony/snn1/applied.qs"     
      [12] "test_data_dir/test_user/default/single-cell/mock_10x_integrated_harmony/snn1/clusters.qs"    
      [13] "test_data_dir/test_user/default/single-cell/mock_10x_integrated_harmony/snn1/scseq_sample.qs"
      [14] "test_data_dir/test_user/default/single-cell/mock_10x_integrated_harmony/snn1/summed.qs"      
      [15] "test_data_dir/test_user/default/single-cell/mock_10x_integrated_harmony/snn_graph.qs"        
      [16] "test_data_dir/test_user/default/single-cell/mock_10x_integrated_harmony/species.qs"          
      [17] "test_data_dir/test_user/default/single-cell/mock_10x_integrated_harmony/tlogs.tenx"          

---

    Code
      integrated_clusters
    Output
      [1] "1" "2" "3" "4" "5"

