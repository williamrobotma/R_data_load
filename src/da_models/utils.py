"""Utility functions for DA models."""

import os

import h5py
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import paired_cosine_distances


def get_metric_ctp(metric_name):
    """Get metric class.

    Args:
        metric_name (str): Metric name.

    Returns:
        Callable Metric.

    """

    if metric_name == "cos":

        def avg_cosine_distance(x, y):
            return paired_cosine_distances(x, y).mean()

        return avg_cosine_distance

    raise ValueError(f"Unknown metric: {metric_name}")


def dict_to_lib_config(args_dict):
    return None


class ModelWrapper:
    @classmethod
    def configurate(cls, lib_config):
        pass

    def __init__(self, model_or_model_dir, name="final_model"):
        prediction_dir = os.path.join(model_or_model_dir, "final_model")
        with h5py.File(os.path.join(prediction_dir, "st_predictions.h5"), "r") as f:
            self.predictions = f["predictions"][()]

        col_names = pd.read_csv(
            os.path.join(prediction_dir, "pred_columns.csv"), header=None, index_col=None
        )
        row_names = pd.read_csv(
            os.path.join(prediction_dir, "pred_rows.csv"), header=None, index_col=None
        )

        self.predictions = pd.DataFrame(self.predictions, index=row_names, columns=col_names)

    def get_predictions(self, input, source_encoder=False):
        return self.model.predict(input)

    def get_embeddings(self, input, source_encoder=False):
        return self.embs.predict(input)


def get_best_params_file(model_name, dset, sc_id, st_id, configs_dir="configs"):
    return
    # pattern = os.path.join(
    #     "model",
    #     model_name,
    #     dset,
    #     f"{sc_id}_{st_id}",
    #     "**",
    #     "reverse_val_best_epoch.csv",
    # )

    # results = []
    # for rv_result_path in glob.glob(pattern, recursive=True):
    #     results.append(pd.read_csv(rv_result_path, index_col=0))

    # results_df = pd.concat(results, axis=0)
    # best_hp = results_df["val"].idxmin()
    # config_fname = results_df.loc[best_hp, "config_fname"]
    # with open(os.path.join(configs_dir, model_name, config_fname), "r") as f:
    #     config = yaml.safe_load(f)

    # lib_params = config["lib_params"]
    # data_params = config["data_params"]
    # model_params = config["model_params"]

    # torch_seed = lib_params.get("manual_seed")
    # lib_seed_path = str(torch_seed) if "manual_seed" in lib_params else "random"

    # model_folder = data_loading.get_model_rel_path(
    #     model_name,
    #     model_params["model_version"],
    #     lib_seed_path=lib_seed_path,
    #     **data_params,
    # )

    # model_folder

    # model_path = os.path.join(
    #     "model",
    #     model_folder,
    #     "advtrain",
    #     "samp_split" if data_params.get("samp_split") else "",
    #     "final_model.pth",
    # )
    # checkpoint = torch.load(model_path)

    # try:
    #     epoch = checkpoint["epoch"]
    # except KeyError:
    #     epoch = checkpoint.get("iters")

    # if int(epoch) != int(results_df.loc[best_hp, "best_epoch"]):
    #     raise ValueError("Epoch mismatch")

    # return config_fname, results_df, best_hp
