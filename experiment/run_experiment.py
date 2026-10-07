#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
run_experiment.py

Runs one session of the cross-modal transitivity experiment (loudness,
brightness and vibration strength). Participants adjust the intensity of a
target stimulus until it matches the standard, or is p times as intense
(production factor p).

Start from the experiment folder with
    python run_experiment.py
and enter the participant alias (e.g. 01) and session number (01, 02, ...).

Session flow:
  1. Demographics (session 01 only, not for alias "test")
  2. 30 s dark adaptation, tutorial, key overview
  3. Training trials, then all experimental blocks with pauses in between
  4. Session file for the next session is created by calling
     session/create_next_session_file_call_from_python.R

Input:  session/subj_<alias>/subj_<alias>session_<ses>.txt  (trial list)
Output: data/subj_<alias>/subj<alias>ses<ses>_<timestamp>.txt  (trial data)
        session/subj_<alias>/standard_match.txt  (matches of all sessions,
                                                  used for start values)
        data/logFile_subj_<alias>.txt  (log of problematic trials)
        data/demographics.txt

Created:  2024/11/26, AO
Last mod: 2024/12/21, AO; 2025/10/21, GM; 2026/09/30, Ina Kohler
"""
import time


if __name__ == "__main__":
    # ── Participant info ─────────────────────────────────────────────────────
    expinfo = dict(alias = input("Alias (z.B. 01): "),
                   session = input("Session (01, 02, ...): "),
                   tstamp = time.strftime("%Y%m%d_%H%M", time.localtime()))
    # Demographics separately (and only in session 1)
    if expinfo["session"] == "01" and expinfo["alias"] != "test":
        demographics = dict(gender = input("Gender (f/m/d): "),
                            age = input("Alter: "))
        with open("data/demographics.txt", "a") as demoFile:
            demoFile.write(demographics['gender'] + '\t' + str(demographics['age']) + '\n')

    # PsychoPy and the tactor connection are only loaded after the console
    # input
    import pandas as pd
    from psychopy import core, visual, event
    from trialFunction import trial
    from global_setup import globalSetup
    import subprocess
    import os
    # Imported only to connect to the tactor controller
    # (see vib_exp_v2_successive.py)
    import vib_exp_v2_successive

    import json
    from datetime import timedelta
    from pprint import pprint

    # Window, stimuli and experimental parameters
    specs = globalSetup()

    # ── Files ────────────────────────────────────────────────────────────────
    # Session file created by createSession.R. Columns (itertuples index):
    #   [1] id  [2] block  [3] trial_type  [4] standard_modality  [5] standard
    #   [6] target_modality  [7] start_target  [8] prod_factor
    #   [9] presentation_order
    with open("session/subj_" + expinfo["alias"] + "/subj_" + expinfo["alias"] + "session_" + expinfo["session"] + ".txt") as sesFile:
        stimuli = pd.read_table(sesFile, header=0, sep= " ")
    os.makedirs("data/subj_" + expinfo["alias"], exist_ok=True)
    dataFile = open("data/subj_" + expinfo["alias"] + "/subj" + expinfo["alias"] + "ses" + expinfo["session"] + "_" + time.strftime("%Y%m%d-%H%M%S") + ".txt", "w")
    dataFile.write("#" + str(expinfo) + "\n")
    dataFile.write("id block trial_type target_modality standard start match p keys rt exptime lightRGB standard_modality\n")
    # Matches of all sessions, read by createSession.R to set the target
    # start values of the next session
    stanMatch = open("session/subj_" + expinfo["alias"] + "/standard_match.txt", "a")
    if expinfo["session"] == "01":
        stanMatch.write("session standard target_modality standard_modality prod_factor match\n")
    # Open file to log problematic trials
    logFile = open("data/logFile_subj_" + expinfo["alias"] + ".txt", "a")
    logFile.write("targetModality standard compIntensity p\n")

    # Block length and number of blocks are taken from the session file
    specs.blockLength = len(stimuli[stimuli["block"].astype(str) == "1"].index)
    specs.nblocks = len(set(stimuli[stimuli["block"].astype(str) != "training"]["block"]))

    # ── Adaptation and tutorial ──────────────────────────────────────────────
    specs.present_text("Willkommen! Mit 'Enter' geht es los.")
    specs.adaptation(30)

    specs.tutorial()
    specs.present_text(specs.key_overview)

    # ── Trials ───────────────────────────────────────────────────────────────
    for stim_tuple in stimuli.itertuples():
        # Short message after training is complete
        if stim_tuple[1] != "training" and specs.trialCount == 0:
            specs.present_text("Training beendet!\n\nMit 'Enter' geht es los!")

        # NOTE: for stim tuples, indices start at 1 (index 0 is the row index)
        if stim_tuple[3] == "basic":
            if stim_tuple[1] == "training":
                # Show the full instruction for the production factor p
                # before each training trial
                specs.win.flip()
                core.wait(0.2)
                visual.TextStim(specs.win, text=specs.p_instruction(stim_tuple[6],
                                                                    stim_tuple[4],
                                                                    stim_tuple[8],
                                                                    stim_tuple[9])[0],
                                height=20, colorSpace="rgb255", color=60).draw()
                specs.win.flip()
                core.wait(2)
            standard, start, match, keys, rt, standard_modality = trial(specs,
                        stim_tuple[6], stim_tuple[5], stim_tuple[7],
                        stim_tuple[8], stim_tuple[9], stim_tuple[4],
                        logFile)
            if stim_tuple[1] == "training":
                visual.TextStim(specs.win, "Finale Einstellung gespeichert", height=20,
                                colorSpace="rgb255", color=60).draw()
                specs.win.flip()
                core.wait(0.5)
            # Store the match under the id of the successive trial, which
            # uses it as its standard (<block>.1.<n> -> <block>.2.<n>)
            specs.basicDict[specs.find_successive_id(stim_tuple[1])] = match
        elif stim_tuple[3] == "successive":
            # Standard is the match from the corresponding basic trial
            standard, start, match, keys, rt, standard_modality = trial(specs,
                        stim_tuple[6], specs.basicDict[stim_tuple[1]],
                        stim_tuple[7], stim_tuple[8],
                        stim_tuple[9], stim_tuple[4], logFile)
        # Save data
        dataFile.write(f"{stim_tuple[1]} {stim_tuple[2]} {stim_tuple[3]} {stim_tuple[6]} {standard} {start} {match} {stim_tuple[8]} {keys} {round(rt, 2)} {round(specs.total_timer.getTime(), 2)} {specs.light_stim.fillColor[0]} {standard_modality}\n")
        # Experimental trials: save match for the next session file and
        # pause after each block
        if stim_tuple[1] != "training":
            stanMatch.write(f"{expinfo['session']} {standard} {stim_tuple[6]} {stim_tuple[4]} {stim_tuple[8]} {match}\n")
            specs.trialCount += 1
            if specs.trialCount % specs.blockLength == 0:
                specs.block_pause(int(specs.trialCount / specs.blockLength), specs.nblocks)

    dataFile.close()
    stanMatch.close()
    logFile.close()

    # ── End of session ───────────────────────────────────────────────────────
    specs.present_text("Diese Sitzung ist beendet.\nVielen Dank f\u00fcr deine Teilnahme!")
    specs.win.close()
    # Create session file for the next session. Path to Rscript.exe must
    # match the R installation on the lab computer.
    subprocess.run([
        "C:\Program Files\R\R-4.4.3\\bin\Rscript.exe",
        "session\create_next_session_file_call_from_python.R",
        expinfo["alias"],          # passed as args[1] in R
        str(int(expinfo["session"]) + 1)   # next session number
    ], shell=True)
    core.quit()
