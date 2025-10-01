#%%
import dask.dataframe as dd
# import dask
from dask_jobqueue import SLURMCluster
from dask.distributed import Client
from pathlib import Path
import shutil

class DaskController:
    def __init__(
        self,
        local: bool = False,
        n_workers: int = 10,
        n_workers_adapt_max: int = 10,
        queue: str = 'short',
        memory: str = '4GB',
        walltime: str = '01:00:00',  # Format 'HH:MM:SS'
        local_directory: str | Path = Path('/home/lag36/scratch/lag36/dask'),
    ):
        self.local = local
        self.n_workers = n_workers
        self.n_workers_adapt_max = n_workers_adapt_max
        self.queue = queue
        self.memory = memory
        self.walltime = walltime
        self.local_directory = Path(local_directory)
        self.client = None
        self.cluster = None
    
    def start_dask(self):
        if self.local:
            self.client = Client()
        else:
            self.cluster = SLURMCluster(
                n_workers=self.n_workers,
                queue=self.queue,          # Replace with your SLURM queue/partition
                death_timeout='60',
                cores=1,
                processes=1,
                memory=self.memory,
                walltime=self.walltime,
                local_directory=self.local_directory,
                log_directory=self.local_directory,
            )
            self.cluster.adapt(minimum=self.n_workers, maximum=self.n_workers_adapt_max)
            self.client = Client(self.cluster)
        print(f"Dask Dashboard is running at: {self.client.dashboard_link}")
        
    def shutdown(self):
        if self.client:
            self.client.shutdown()
        if self.cluster:
            self.cluster.close()
        if self.local_directory.exists():
            shutil.rmtree(self.local_directory)            

#%%
dask_controller = DaskController(
    local=False,
    n_workers=20,
    n_workers_adapt_max=30,
    queue='short',
    memory='8GB',
    walltime='03:00:00',
    local_directory='/home/lag36/scratch/lag36/dask',
)

# #%%
# dask_controller.start_dask()
# #%%
# dask_controller.shutdown()
# #%%
# import submarlin_postprocessing.filepaths as filepaths
# df_barcodes = dd.read_parquet(filepaths.final_barcodes_df_merged_filenames['lLAG08'], engine='pyarrow')


#%%
# df_barcodes.columns
# # %%
# df_barcodes.reset_index().drop_duplicates(subset='Multi-Experiment Phenotype Trenchid').head()