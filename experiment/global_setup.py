"""
global_setup.py

Class globalSetup holding the window, stimuli and experimental parameters,
and all functions for presenting stimuli, texts, pauses and the tutorial.

Stimuli:
  visual   : grey circle on black background (light_stim; light_stim_2 is
             the second circle in visual -> visual trials)
  auditory : 1000 Hz sine tone
  tactile  : 500 ms vibration pulse via the tactor controller
             (TactorInterface.dll)
Background: brown noise (smoothed_brown_noise_8hrs.mp3) plays during the
whole session to mask the sound of the tactor.

All intensities are in device units (see trialFunction.py).

Last mod: 2024/12/21, AO; 2025/03/17, KN; 2025/10/21, GM;
          2026/09/30, Ina Kohler
"""

from psychopy import visual, core, event
from psychopy import prefs
prefs.hardware['audioLib'] = ['PTB']
prefs.hardware['audioDevice'] = [1.0]
from psychopy import sound
import numpy as np
import pandas as pd
import datetime
import vib_exp_v2_successive
from conversions import visualscale_to_rgb
from ctypes import *
import pygame
# Library for the tactor controller (connected in vib_exp_v2_successive)
lib = CDLL('TactorInterface.dll')

class globalSetup:

    def __init__(self):
        # ── Window ───────────────────────────────────────────────────────────
        self.win = visual.Window([2048,1536], units='pix', color=0, screen=1,
                            colorSpace="rgb255", allowGUI=False, fullscr=True)
        self.fixation = visual.TextStim(self.win, text='+', units='pix',
                            pos=(0,0), height=40, colorSpace='rgb255', color=60)

        # ── Auditory stimulus ────────────────────────────────────────────────
        # Tone amplitude = 10**((intensity - offset_spl)/20) * sin_wave
        self.offset_spl = 94
        self.stimulus_duration = 0.50
        # Loaded but not used; the tone is a 1000 Hz sine
        self.pink_noise = np.loadtxt("stimuli/pinkNoiseArray.txt")
        t = np.linspace(0, 1.0, int(44100*1.0), endpoint=False)
        self.sin_wave = 0.5*np.sin(2*np.pi*1000*t)

        self.sound_stim = sound.Sound(value=self.sin_wave, sampleRate=44100,
                            secs=self.stimulus_duration, stereo=True)

        # ── Visual stimuli ───────────────────────────────────────────────────
        self.light_stim = visual.Circle(self.win, color=70, pos=(0, 0),
                            size=121, edges=128, colorSpace="rgb255")
        self.light_stim_2 = visual.Circle(self.win, color=70, pos=(0, 0),
                            size=121, edges=128, colorSpace="rgb255")

        # ── Response keys ────────────────────────────────────────────────────
        # Step sizes in device units. Lower row (a, s, d) decreases, upper
        # row (w, e, r) increases the target; step size grows from left to
        # right.
        self.auditory_steps = dict(a=-1, s=-4, d=-6, w=1, e=4, r=6)
        self.visual_steps = dict(a=-1, s=-3, d=-6, w=1, e=3, r=6)
        self.vibrational_steps = dict(a=-1, s=-5, d=-20, w=1, e=5, r=20)
        self.allowed_keys = ['q', 'space', 'a', 's', 'd', 'w', 'e', 'r']
        self.key_overview = """Tastenübersicht:
                    \n\nObere Reihe:\nw (+), e (++), r (+++)
                    \n\nUntere Reihe:\na (-), s (--), d (---)
                    \nWenn du zufrieden mit deiner Einstellung bist, drücke die Leertaste.
                    \n\nWeiter mit 'Enter'."""
        # Short instruction shown above the stimuli during each trial
        self.p_instr = visual.TextStim(self.win, "", height=20,
                            colorSpace="rgb255", pos=(0, 170), color=60)

        # ── Timing (s) ───────────────────────────────────────────────────────
        self.fixation_time = 0.5
        self.isi = 0.05

        # ── Session state ────────────────────────────────────────────────────
        # blockLength is overwritten in run_experiment.py from the session
        # file
        self.blockLength = 80
        self.rt_timer = core.Clock()
        self.total_timer = core.Clock()
        # Matches of basic trials, keyed by the id of the successive trial
        # that uses them as standard
        self.basicDict = dict()
        self.trialCount = 0

        # ── Background noise ─────────────────────────────────────────────────
        pygame.mixer.init()
        pygame.mixer.music.load('smoothed_brown_noise_8hrs.mp3')
        pygame.mixer.music.set_volume(0.1)  # 0.0 to 1.0
        pygame.mixer.music.play(-1)
        print("init aus global setup beendet")


    # ── Presentation ─────────────────────────────────────────────────────────

    def present_fixation(self):
        """Present fixation cross for 300ms"""
        self.fixation.draw()
        self.win.flip()
        core.wait(0.3)
        self.win.flip()

    def present_text(self, text):
        """Create a TextStim object with the input text and present it on the
        current window until return key is pressed."""
        visual.TextStim(self.win, text, height=20, colorSpace="rgb255", color=60,
                        alignText="left").draw()
        self.win.flip()
        event.waitKeys(keyList='return')

    def present_sound(self):
        """Play sound stimulus for a duration of global variable 'stimulus
        duration'"""
        snd_timer = core.CountdownTimer(self.stimulus_duration)
        while snd_timer.getTime() > 0:
            self.sound_stim.play()
            core.wait(min(self.stimulus_duration, snd_timer.getTime()))
            self.sound_stim.stop()

    def present_light(self):
        """Present light stimulus for a duration of global variable 'stimulus
        duration'"""
        self.light_stim.draw()
        self.win.flip()
        core.wait(self.stimulus_duration)
        self.win.flip()

    def present_vibration(self):
        """Start a 500 ms vibration pulse at intensity self.vib_intensity"""
        lib.ChangeGain(c_int(0), c_int(1), c_int(self.vib_intensity), c_int(0))
        lib.Pulse(c_int(0), c_int(1), c_int(500), c_int(0))


    def present_stimuli(self, target_modality, standard_modality, presentation_order, new_intensity):
        """Present standard and target consecutively after a fixation cross.

        presentation_order 1 = target first, 2 = standard first. Each branch
        below covers one modality pair in both orders. The first stimulus is
        presented for stimulus_duration, followed by the ISI. Keys can be
        pressed from 200 ms after onset of the second stimulus; a key press
        ends the presentation early.

        new_intensity is the current target intensity, only needed for the
        vibration (other targets are updated in adjust_intensity()).

        Returns the list of keys pressed during the second stimulus, or None.
        """
        # Fixation
        self.fixation.draw()
        self.win.flip()
        core.wait(self.fixation_time)
        self.win.flip()

        # Stimuli
        # light -> sound
        if ((target_modality == "visual" and standard_modality == "auditory" and presentation_order == 1) or
            (target_modality == "auditory" and standard_modality == "visual" and presentation_order == 2)):
            self.present_light()
            core.wait(self.isi)
            self.sound_stim.play()
            core.wait(0.2)
            # allow response after 200ms
            keylist = event.waitKeys(maxWait=self.stimulus_duration-0.2,
                                     keyList=self.allowed_keys)
            self.sound_stim.stop()
            self.win.flip()
        # sound -> light
        elif ((target_modality == "auditory" and standard_modality == "visual" and presentation_order == 1) or
             (target_modality == "visual" and standard_modality == "auditory" and presentation_order == 2)):
            self.present_sound()
            core.wait(self.isi)
            self.light_stim.draw()
            self.win.flip()
            core.wait(0.2)
            # allow response after 200ms
            keylist = event.waitKeys(maxWait=self.stimulus_duration-0.2,
                                     keyList=self.allowed_keys)
            self.win.flip()

        # vibration -> light
        elif ((target_modality == "tactile" and standard_modality == "visual" and presentation_order == 1) or
             (target_modality == "visual" and standard_modality == "tactile" and presentation_order == 2)):
            self.vib_intensity = new_intensity if target_modality == "tactile" else self.vib_standard
            self.present_vibration()
            core.wait(self.isi)
            self.light_stim.draw()
            self.win.flip()
            core.wait(0.2)
            # allow response after 200ms
            keylist = event.waitKeys(maxWait=self.stimulus_duration-0.2,
                                     keyList=self.allowed_keys)
            self.win.flip()

        # light -> vibration
        elif ((target_modality == "visual" and standard_modality == "tactile" and presentation_order == 1) or
             (target_modality == "tactile" and standard_modality == "visual" and presentation_order == 2)):
            self.present_light()
            core.wait(self.isi)
            self.vib_intensity = new_intensity if target_modality == "tactile" else self.vib_standard
            self.present_vibration()
            core.wait(0.2)
            # allow response after 200ms
            keylist = event.waitKeys(maxWait=self.stimulus_duration-0.2,
                                     keyList=self.allowed_keys)
            self.win.flip()

        # sound -> vibration
        elif ((target_modality == "auditory" and standard_modality == "tactile" and presentation_order == 1) or
             (target_modality == "tactile" and standard_modality == "auditory" and presentation_order == 2)):
            self.present_sound()
            core.wait(self.isi)
            self.vib_intensity = new_intensity if target_modality == "tactile" else self.vib_standard
            self.present_vibration()
            core.wait(0.2)
            # allow response after 200ms
            keylist = event.waitKeys(maxWait=self.stimulus_duration-0.2,
                                     keyList=self.allowed_keys)
            self.win.flip()

        # vibration -> sound
        elif ((target_modality == "tactile" and standard_modality == "auditory" and presentation_order == 1) or
             (target_modality == "auditory" and standard_modality == "tactile" and presentation_order == 2)):
            self.vib_intensity = new_intensity if target_modality == "tactile" else self.vib_standard
            self.present_vibration()
            core.wait(self.isi)
            self.sound_stim.play()
            core.wait(0.2)
            # allow response after 200ms
            keylist = event.waitKeys(maxWait=self.stimulus_duration-0.2,
                                     keyList=self.allowed_keys)
            self.sound_stim.stop()
            self.win.flip()

        # light -> light, target (light_stim) first
        elif ((target_modality == "visual" and standard_modality == "visual" and presentation_order == 1)):
            self.light_stim.draw()
            self.win.flip()
            core.wait(self.stimulus_duration)
            self.win.flip()
            core.wait(self.isi)
            self.light_stim_2.draw()
            self.win.flip()
            core.wait(0.2)
            # allow response after 200ms
            keylist = event.waitKeys(maxWait=self.stimulus_duration-0.2,
                                     keyList=self.allowed_keys)
            self.win.flip()
        # light -> light, standard (light_stim_2) first
        elif ((target_modality == "visual" and standard_modality == "visual" and presentation_order == 2)):
            self.light_stim_2.draw()
            self.win.flip()
            core.wait(self.stimulus_duration)
            self.win.flip()
            core.wait(self.isi)
            self.light_stim.draw()
            self.win.flip()
            core.wait(0.2)
            # allow response after 200ms
            keylist = event.waitKeys(maxWait=self.stimulus_duration-0.2,
                                     keyList=self.allowed_keys)
            self.win.flip()

        self.fixation.draw()
        self.win.flip()
        return keylist

    def adjust_intensity(self, target_modality, new_intensity):
        """Clamp the target intensity to the modality's limits and update the
        target stimulus. Returns the clamped intensity.

        Limits (device units): auditory 15-100, visual 1-49, tactile 10-255.
        """
        if target_modality == "auditory":
            if new_intensity < 15:
                new_intensity = 15
            elif new_intensity > 100:
                new_intensity = 100
            self.sound_stim.setSound(value=10**((new_intensity-self.offset_spl)/20) *
                    self.sin_wave)  # previously self.pink_noise
        elif target_modality == "visual":
            if new_intensity < 1:
                new_intensity = 1
            elif new_intensity > 49:
                new_intensity = 49
            self.light_stim.lineColor = visualscale_to_rgb(new_intensity)
            self.light_stim.fillColor = visualscale_to_rgb(new_intensity)
        elif target_modality == "tactile":
            if new_intensity < 10:
                new_intensity = 10
            elif new_intensity > 255:
                new_intensity = 255
            self.vib_intensity = new_intensity

        return new_intensity

    def find_successive_id(self, basicID):
        """For basic trial IDs, return ID of respective successive trial"""
        return basicID.replace(".1.", ".2.")

    # ── Instructions and pauses ──────────────────────────────────────────────

    def p_instruction(self,
                      target_modality,
                      standard_modality,
                      p,
                      presentation_order):
        """Create the instruction for production factor p.

        Returns a tuple: the full instruction sentence (shown before training
        trials) and the short name of the target (shown during each trial,
        e.g. "Ton"; "1. Kreis" / "2. Kreis" for visual -> visual).
        """
        if p == 1:
            prod = "genauso"
        else:
            prod = str(p) + "-mal so"
        if target_modality == "visual":
            if standard_modality == "visual":
                iw = ("Helligkeit", "Kreises", prod, "hell", str(presentation_order) + ". Kreis")
            else:
                iw = ("Helligkeit", "Kreises", prod, "hell", "Kreis")
        elif target_modality == "auditory":
            iw = ("Lautstärke", "Ton", prod, "laut", "Ton")
        elif target_modality == "tactile":
            iw = ("Stärke", "Pulses", prod, "stark", "Puls")

        if standard_modality == "visual":
            iw = iw + ("Helligkeit", "Kreises")
        elif standard_modality == "auditory":
            iw = iw + ("Lautstärke", "Ton")
        elif standard_modality == "tactile":
            iw = iw + ("Stärke", "Pulses")
        return (f"""Mache die {iw[0]} des {iw[1]} {iw[2]} {iw[3]} wie die {iw[5]} des {iw[6]}.""", iw[4])

    def block_pause(self, thisblock, nblocks):
        """Pause screen after each block except the last"""
        if int(thisblock) == nblocks:
            return
        else:
            self.present_text("Pause! " + str(thisblock) + " von " + str(nblocks) + " Blöcken geschafft!\nWeiter geht's mit der Enter-Taste wenn du bereit bist.")
            self.present_text(self.key_overview)

    def adaptation(self, duration=300):
        """Adaptation to darkness with countdown (duration in s)"""
        adaptation_timer = core.CountdownTimer(duration)
        while adaptation_timer.getTime() > 0:
            visual.TextStim(self.win, "Bevor du beginnst sollen sich deine Augen zuerst an die Dunkelheit gewöhnen. Mache deine Augen also kurz zu und zähle langsam bis 30...:\n\n " + str(datetime.timedelta(seconds=round(adaptation_timer.getTime(), 0))),
                height=20, colorSpace="rgb255", color=60, alignText="left").draw()
            self.win.flip()
        self.present_text("Adaptation beendet. Mit 'Enter' geht es los!")

    # ── Tutorial ─────────────────────────────────────────────────────────────
    # test_* functions let the participant adjust one stimulus freely until
    # space is pressed.

    def test_sound(self, new_intensity):
        """Demo of sound stimulus adjustment"""
        new_intensity = self.adjust_intensity("auditory", new_intensity)
        self.win.flip()
        self.present_sound()
        self.win.flip()
        key = event.waitKeys(keyList=self.allowed_keys).pop()
        if key == 'space':
            self.win.flip()
            return
        else:
            new_intensity += self.auditory_steps[key]
            self.test_sound(new_intensity)

    def test_light(self, new_intensity):
        """Demo of light stimulus adjustment. NOTE: uses the auditory step
        sizes, which differ slightly from the visual ones."""
        new_intensity = self.adjust_intensity("visual", new_intensity)
        self.present_light()
        self.win.flip()
        key = event.waitKeys(keyList=self.allowed_keys).pop()
        if key == 'space':
            self.win.flip()
            return
        else:
            new_intensity += self.auditory_steps[key]
            self.test_light(new_intensity)

    def test_vibration(self, new_intensity):
        """Demo of vibration stimulus adjustment"""
        new_intensity = self.adjust_intensity("tactile", new_intensity)
        self.vib_intensity = new_intensity
        self.present_vibration()
        key = event.waitKeys(keyList=self.allowed_keys).pop()
        if key == 'space':
            self.win.flip()
            return
        else:
            new_intensity += self.vibrational_steps[key]
            self.test_vibration(new_intensity)


    def tutorial(self):
        """Instructions and tutorial for experimental procedure: example of
        each stimulus, then free adjustment of sound, light and vibration"""
        # Show stimuli
        self.present_text("Los geht es mit einer Instruktion.\n\nÜber die Kopfhörer wird in jedem Durchgang ein Geräusch gespielt oder auf dem Bildschirm wird ein grauer Kreis angezeigt oder du spuerst eine Vibration am Unterarm.\n\nWeiter geht es jetzt mit 'Enter', dann hörst du ein Beispiel-Geräusch, siehst einen Beispiel-Kreis und spÜrst eine Beispiel-Vibration.")
        self.win.flip()
        core.wait(0.3)
        self.sound_stim.play()
        core.wait(self.stimulus_duration)
        self.sound_stim.stop()
        self.win.flip()
        core.wait(self.isi)
        self.present_light()
        core.wait(self.isi)
        self.vibrationstest()
        core.wait(0.5)
        # Explain adjustment of stimuli
        self.present_text("Deine Aufgabe in jedem Durchgang wird sein, entweder die Lautstärke des Geräuschs, die Helligkeit des Kreises oder die Stärke der Vibration anzupassen.\n\nDies kannst du mit den Tasten der Tastatur vor dir machen.\n\nWeiter mit 'Enter'.")
        self.present_text("Probiere zuerst aus, das Geräusch zu verändern.\n\nDie Tasten der oberen Reihe machen es lauter und die der unteren Reihe machen es leiser.\n\nWeiter mit 'Enter'.")
        self.present_text('''Verändere jetzt die Lautstärke des Tons.\n\nGroße Schritte (--- oder +++) ➞ Zeigefinger \nmittlere Schritte (-- oder ++) ➞ Mittelfinger \nkleine Schritte (- oder +) ➞ Ringfinger.
                          \n\nWenn du zufrieden mit deiner Einstellung bist, drücke die Leertaste, um den Durchgang zu beenden.\n\nMit 'Enter' geht es weiter.''')
        self.test_sound(50)
        self.present_text("""Probiere nun das gleiche mit dem Licht.\nGroße Schritte (--- oder +++) ➞ Zeigefinger \nmittlere Schritte (-- oder ++) ➞ Mittelfinger \nkleine Schritte (- oder +) ➞ Ringfinger.
                \n\nWenn du zufrieden mit deiner Einstellung bist, drücke die Leertaste, um den Durchgang zu beenden.\n\nMit 'Enter' geht es weiter.""")
        self.test_light(60)
        self.present_text("""Probiere nun das gleiche mit der Vibration.\nGroße Schritte (--- oder +++) ➞ Zeigefinger \nmittlere Schritte (-- oder ++) ➞ Mittelfinger \nkleine Schritte (- oder +) ➞ Ringfinger
                \n\nWenn du zufrieden mit deiner Einstellung bist, drücke die Leertaste, um den Durchgang zu beenden.\n\nMit 'Enter' geht es weiter.""")
        self.test_vibration(180)
        self.present_text("Sehr gut!\nIn jedem Durchgang bekommst du eine Anweisung, welchen der 3 Reize - Ton, Licht oder Vibration - du einstellen sollst, und auch wie stark.\n\n[Enter]")
        self.present_text("Falls du jetzt noch Fragen hast, wende dich bitte an die Versuchsleitung.\n\nWenn du 'Enter' drückst, geht es los mit ein paar Trainingsdurchgängen.")


    def vibrationstest(self):
        """Example vibration in the tutorial; also checks that the tactor
        is working"""
        dev = c_int(0)
        tact1 = c_int(1)
        # Base settings
        current_intensity = 180  # Starting intensity
        duration = 500  # Duration of each vibration in ms
        lib.ChangeGain(dev, tact1, c_int(current_intensity), c_int(0))
        lib.Pulse(dev, tact1, c_int(duration), c_int(0))
        print("Testfunktion durchlaufen")
