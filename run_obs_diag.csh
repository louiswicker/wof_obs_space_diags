#!/bin/csh -vx
#==================================================================

source ~/.tcshrc

# 09 May
#python obs_seq_collate.py -d "/scratch/wof/realtime/20170509/201705*" -f "obs_seq.final*" -p 09May_OAN17 > log&
#python obs_seq_collate.py -d "/scratch/wof/realtime/20170509/EXPS/RLT/DART_OUTPUT" -f "obs_seq.final*" -p 09May_RLT > log& 


# 16 May
python obs_seq_collate.py -d "/scratch/wof/realtime/20170516/201705*" -f "obs_seq.final*" -p 16May_OAN18  >& 16May.log &
#python obs_seq_collate.py -d "/scratch/wof/realtime/20170516/EXPS/RLT/DART_OUTPUT" -f "obs_seq.final*" -p 16May_RLT > log& 

# 17 May
#python obs_seq_collate.py -d "/scratch/wof/realtime/20170517/201705*" -f "obs_seq.final*" -p 17May_OAN17  > log&
#python obs_seq_collate.py -d "/scratch/wof/realtime/20170517/EXPS/RLT/DART_OUTPUT" -f "obs_seq.final*" -p 17May_RLT  > log&

# 18 May
python obs_seq_collate.py -d "/scratch/wof/realtime/20170518/201705*" -f "obs_seq.final*" -p 18May_OAN18 >& 18May.log & 
#python obs_seq_collate.py -d "/scratch/wof/realtime/20170518/EXPS/RLT/DART_OUTPUT" -f "obs_seq.final*" -p 18May_RLT 

# 23 May
#python obs_seq_collate.py -d "/scratch/wof/realtime/20170523/201705*" -f "obs_seq.final*" -p 23May_OAN17 > log&
#python obs_seq_collate.py -d "/scratch/wof/realtime/20170523/EXPS/RLT/DART_OUTPUT" -f "obs_seq.final*" -p 23May_RLT > log& 

# 27 May
#python obs_seq_collate.py -d "/scratch/wof/realtime/20170527/201705*" -f "obs_seq.final*" -p 27May_OAN17  > log&
#python obs_seq_collate.py -d "/scratch/wof/realtime/20170527/EXPS/RLT/DART_OUTPUT" -f "obs_seq.final*" -p 27May_RLT  > log&

