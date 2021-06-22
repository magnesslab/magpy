import magpy.settings as settings
import magpy.pipeline as pipeline
# import magpy.velocity as velocity
import magpy.bulkseq as bulkseq

from magpy.functions import *
from magpy.gene_set_enrichment import *
from magpy.io_functions import *
from magpy.evaluate_hashtags import *
from magpy.diffx_testing import de_test


# Added 5-Jun-2021 (wolfram)
# If running in ipython, add additional garbage collection following exectution of any cell
# This fixes an Ipython/AnnData memory leak bug where adata objects are not removed from memory 
try: 
    ip = get_ipython()
    is_collecting = False
    for function in ip.events.__dict__['callbacks']['post_run_cell']:
        if function.__name__ == '_gccollect': is_collecting = True
    if not is_collecting: 
        import gc
        import psutil
        def _gccollect(): 
#        	if psutil.virtual_memory()[2] > 50:
#        		print("Server RAM usage exceeded 50\% after cell execution, please notify an admin.")
#        		print("If repeatedly overwring adata variables, consider using gc.collect().")
        	gc.collect()
        ip.events.register('post_run_cell',_gccollect)
        
except: print("Ipython not running, skipping gc fix.")
