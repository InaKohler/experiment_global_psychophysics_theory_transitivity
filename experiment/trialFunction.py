"""
trialFunction.py

Runs a single matching trial: sets standard and target intensities, presents
the stimulus pair repeatedly and adjusts the target after each key press
until the participant confirms with space.

All intensities are in device units:
  visual   : visualscale, converted to RGB via the lookup table
  auditory : level relative to specs.offset_spl
  tactile  : gain of the tactor controller (0-255)

Created:  2024, AO
Last mod: 2024/12/21, AO; 2025/10/21, GM; 2026/09/30, Ina Kohler
"""

from psychopy import visual, core, sound, event
from conversions import visualscale_to_rgb
from ctypes import *
import vib_exp_v2_successive
# Library for the tactor controller (connected in vib_exp_v2_successive)
lib = CDLL('TactorInterface.dll')

def trial(specs, target_modality, standard_value, start_target, p, pres_order,
          standard_modality, logFile):
    """Run one matching trial.

    specs             : globalSetup instance (window, stimuli, parameters)
    target_modality   : "visual", "auditory" or "tactile"
    standard_value    : intensity of the standard
    start_target      : starting intensity of the target
    p                 : production factor (target = p times the standard)
    pres_order        : 1 = target presented first, 2 = standard first
    standard_modality : "visual", "auditory" or "tactile"
    logFile           : open file for logging problematic trials

    Returns standard_value, start_target, final target intensity (match),
    string of all pressed keys, response time in s and standard_modality.
    """
    # ── Set stimulus intensities ─────────────────────────────────────────────
    # Target first, then standard. For visual -> visual trials the standard
    # is shown on the second circle (light_stim_2).
    if target_modality == "visual":
        specs.light_stim.lineColor = visualscale_to_rgb(start_target)
        specs.light_stim.fillColor = visualscale_to_rgb(start_target)
    elif target_modality == "auditory":
        specs.sound_stim.setSound(10**((start_target-specs.offset_spl)/20) *
                                      specs.sin_wave)
    elif target_modality == "tactile":
        lib.ChangeGain(c_int(0), c_int(1), c_int(start_target), c_int(0))
        specs.vib_intensity = start_target

    if standard_modality == "visual":
        if target_modality == "visual":
            specs.light_stim_2.lineColor = visualscale_to_rgb(standard_value)
            specs.light_stim_2.fillColor = visualscale_to_rgb(standard_value)
        else:
            specs.light_stim.lineColor = visualscale_to_rgb(standard_value)
            specs.light_stim.fillColor = visualscale_to_rgb(standard_value)

    elif standard_modality == "auditory":
        specs.sound_stim.setSound(10**((standard_value-specs.offset_spl)/20) *
                                      specs.sin_wave)
    elif standard_modality == "tactile":
        lib.ChangeGain(c_int(0), c_int(1), c_int(standard_value), c_int(0))
        specs.vib_standard = standard_value

    pressed_keys = ""
    new_intensity = start_target

    # Short instruction shown above the stimuli during the trial, e.g. "Ton x2"
    specs.p_instr.text = specs.p_instruction(target_modality, standard_modality, p, pres_order)[1] + " x" + str(p)
    specs.p_instr.autoDraw = True

    specs.rt_timer.reset()

    # ── Adjustment loop ──────────────────────────────────────────────────────
    # Present the pair, read one key, adjust the target and repeat.
    # Keys pressed while the second stimulus is shown are used directly,
    # otherwise the loop waits for a key.
    while True:
        keylist = specs.present_stimuli(target_modality, standard_modality, pres_order, new_intensity)
        if keylist == None:
            key = event.waitKeys(keyList=specs.allowed_keys).pop()
        else:
            key = keylist.pop()

        # At least one adjustment is required before confirming
        if key == "space" and pressed_keys == "":
            visual.TextStim(specs.win,
                          text="Bitte nimm mindestens eine Einstellung vor!",
                            height=20, colorSpace="rgb255", color=60).draw()
            specs.win.flip()
            core.wait(1)
            specs.win.flip()
            continue

        pressed_keys += key
        if key == "q":
            # Abort experiment
            core.quit()
        elif key == "space":
            # Confirm final adjustment
            response_time = specs.rt_timer.getTime()
            specs.p_instr.autoDraw = False
            return standard_value, start_target, new_intensity, pressed_keys, response_time, standard_modality
        elif key == "h":
            # Log problematic trial. NOTE: "h" is not in specs.allowed_keys,
            # so this branch is never reached.
            logFile.write(f"{target_modality}, {standard_value}, {new_intensity}, {p}\n")
        else:
            # Adjustment keys a/s/d (down) and w/e/r (up), step sizes per
            # modality in globalSetup
            if target_modality == "visual":
                new_intensity += specs.visual_steps[key]
            elif target_modality == "auditory":
                new_intensity += specs.auditory_steps[key]
            elif target_modality == "tactile":
                new_intensity += specs.vibrational_steps[key]
            # Clamp to the modality's limits and update the target stimulus
            new_intensity = specs.adjust_intensity(target_modality, new_intensity)
