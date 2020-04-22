import shutil
import os
import json
import pele_platform.Utilities.Helpers.helpers as hp


def external_adaptive_file(env):
    hp.silentremove([env.ad_ex_temp])
    env.ad_ex_temp = os.path.join(env.pele_dir, os.path.basename(env.adapt_conf))
    shutil.copy(env.adapt_conf, env.ad_ex_temp)
    if env.input:
        with open(env.ad_ex_temp, "r") as json_file:
            data = json.load(json_file)
            data["generalParams"]["initialStructures"] = env.inputs_simulation
        with open(env.ad_ex_temp, "w") as json_file:
            json.dump(data, json_file, sort_keys=True, indent=4)


def external_pele_file(env):
    hp.silentremove([env.pele_exit_temp])
    env.pele_exit_temp = os.path.join(env.pele_dir, os.path.basename(env.confile))
    shutil.copy(env.confile, env.pele_exit_temp)
