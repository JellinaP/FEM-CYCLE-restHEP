"""
Heartbeat and Tone Counting Task Script

A PsychoPy experiment implementing interoceptive (heartbeat counting) and exteroceptive
(tone counting) trials across alternating blocks. Includes a training phase with an
adaptive staircase procedure to calibrate tone volume to each participant's detection
threshold before the main task begins.

Experiment flow:
    1. Training phase - Staircase tone volume calibration (exteroceptive only)
    2. Main task - Alternating blocks of interoceptive and exteroceptive trials,
       with adaptive difficulty adjustments within each condition

Participant groups: NC (natural cycle), OC (oral contraceptive), MA (male)
Each group has condition-specific phases (NC: luteal/follicular, OC: active/placebo, MA: N/A)

Original code by DEO Larsson
Modified by Jellina Prinsen, Dept. of Neurology, MPI CBS (prinsen@cbs.mpg.de)
"""

# ============================================================
# IMPORTS
# ============================================================
import sys
from psychopy import gui, event, visual, core, logging, prefs, sound, parallel
import glob, os, random
from pyglet.window import key
import numpy as np
import time
from datetime import datetime
from ctypes import windll

# ============================================================
# EXPERIMENT PARAMETERS
# ============================================================
ppPresent = 1  # Set to 0 if no parallel port connection (e.g. for testing)

training_running = True      # Set to False to skip training (uses default volume)
n_blocks = range(0, 4)       # Number of blocks; 12 blocks for full task
conditions = (0, 1)          # 0 = interoceptive (heart); 1 = exteroceptive (tones)
n_trials = 3                 # Trials per block; 10 for full task

# Time parameters
intervals = [20, 20, 20, 20, 20, 20, 20, 20, 18, 22]  # 80% at 20s, 20% at +/-2s
iti = 3       # Inter-trial interval (seconds)
final_t = 3   # Post-experiment wait time (seconds)

# ============================================================
# PATHS
# ============================================================
path = os.path.abspath('MINT_Experiment')
rawpath = os.path.abspath('MINT_HCT_BehavioralData')

# ============================================================
# PARTICIPANT GUI
# ============================================================
# First dialog: subjID, Group, Session
subj_inputs = {
    'subjID': '',
    'HormonalGroup': ['NC', 'OC', 'MA'],
    'Session': ['1', '2'],
}
ID_Gui = gui.DlgFromDict(subj_inputs, title="Participant's Information", order=['subjID', 'HormonalGroup', 'Session'])
if not ID_Gui.OK:
    sys.exit('User cancelled. Exiting...')

# Second dialog: Phase (conditional on group)
if subj_inputs['HormonalGroup'] == 'NC':
    phase_inputs = {'Phase': ['luteal', 'follicular']}
    Phase_Gui = gui.DlgFromDict(phase_inputs, title="Phase Selection")
    if not Phase_Gui.OK:
        sys.exit('User cancelled. Exiting...')
    subj_inputs['Phase'] = phase_inputs['Phase']
elif subj_inputs['HormonalGroup'] == 'OC':
    phase_inputs = {'Phase': ['active', 'placebo']}
    Phase_Gui = gui.DlgFromDict(phase_inputs, title="Phase Selection")
    if not Phase_Gui.OK:
        sys.exit('User cancelled. Exiting...')
    subj_inputs['Phase'] = phase_inputs['Phase']
elif subj_inputs['HormonalGroup'] == 'MA':
    subj_inputs['Phase'] = 'NA'

# Get date and timestamp
now = datetime.now()
date_str = now.strftime('%Y-%m-%d')
time_str = now.strftime('%H:%M:%S')

# ============================================================
# PARTICIPANT LIST (group-specific CSV)
# ============================================================
fileName_list = 'Participants_' + str(subj_inputs['HormonalGroup']) + '.csv'

if not os.path.exists(fileName_list):
    part_list = open(fileName_list, 'a')
    part_list.write('subjID,HormonalGroup,Session,Phase,Date,Time\n')
    part_list.write(str(subj_inputs['subjID']) + ',' + str(subj_inputs['HormonalGroup']) + ',' + str(subj_inputs['Session']) + ',' + str(subj_inputs['Phase']) + ',' + date_str + ',' + time_str + '\n')
    part_list.close()
else:
    with open(fileName_list) as f:
        file_content = f.read()
    if str(subj_inputs['subjID']) in file_content and str(subj_inputs['Session']) in file_content:
        sys.exit('WARNING! subjID already exists. Exiting...')
    else:
        part_list = open(fileName_list, 'a')
        part_list.write(str(subj_inputs['subjID']) + ',' + str(subj_inputs['HormonalGroup']) + ',' + str(subj_inputs['Session']) + ',' + str(subj_inputs['Phase']) + ',' + date_str + ',' + time_str + '\n')
        part_list.close()

# ============================================================
# OUTPUT FILE SETUP
# ============================================================
subjID_padded = str(subj_inputs['subjID']).zfill(3)  # Zero-pad to 3 digits
datapath = os.path.join(rawpath, 'sub-' + str(subj_inputs['HormonalGroup']) + subjID_padded, 'sess-' + str(subj_inputs['Session']))

if not os.path.exists(datapath):
    os.makedirs(datapath)

# Build subject file prefix (omit phase from filename for MA group)
subj_prefix = 'sub-' + str(subj_inputs['HormonalGroup']) + subjID_padded
if subj_inputs['HormonalGroup'] != 'MA':
    subj_prefix += '_' + str(subj_inputs['Phase'])

# Interoceptive results
fileName_int_results = subj_prefix + '_task-interoception.csv'
fullName_int_results = os.path.join(datapath, fileName_int_results)
results_int = open(fullName_int_results, 'w')
results_int.write('Trial_Tot,Trial_Interval,Trial_Length,Block,Condition,HB_actual,HB_reported,Accuracy,PAS,BPM\n')

# Exteroceptive results
fileName_ext_results = subj_prefix + '_task-exteroception.csv'
fullName_ext_results = os.path.join(datapath, fileName_ext_results)
results_ext = open(fullName_ext_results, 'w')
results_ext.write('Trial_Tot,Trial_Interval,Trial_Length,Block,Condition,Tone_vol,Tones_Actual,Tones_Reported,Accuracy,PAS,BPM\n')

# Training results
fileName_training_results = subj_prefix + '_task-training.csv'
fullName_training_results = os.path.join(datapath, fileName_training_results)
results_training = open(fullName_training_results, 'w')
results_training.write('Trial_nr,Tones_reported,Tones_actual,Trial_hitRate,Tone_vol,Volume_HighOrLow,PAS\n')

# ============================================================
# PARALLEL PORT AND MARKERS
# ============================================================
if ppPresent:
    from intart_functions import sendParallelTrigger, readParallelTrigger

if ppPresent:
    pport_address_out = 888  # 0x378
    pport_address_in = parallel.ParallelPort(address=0xE010)  # LPT2

    pport = windll.inpoutx64
    pport.Out32(pport_address_out, 0)  # Reset all pins to zero

    # Flip direction bit 5 (pin 7) of control register so input port can read
    pport_address_ctr = parallel.ParallelPort(address=0xE010 + 2)
    if pport_address_ctr.readPin(7) == 0:
        pport_address_ctr.setPin(7, 1)

# Marker codes
mrkrCode_tone = 240
triggerStart = 200
triggerEnd = 205
triggerStartTraining = 100
triggerEndTraining = 105
heartbeatCode = 90

# ============================================================
# WINDOW AND VISUAL STIMULI
# ============================================================
win = visual.Window(
    size=(1100, 1100),
    units='pix',
    color=(0, 0, 0),
    screen=1,
    fullscr=True)

msgpos = (0, 0)
msgheight = 30

# Fixation cross
cross_length = ((0, -3), (0, 3), (0, 0), (-3, 0), (3, 0))
fixation = visual.ShapeStim(win, vertices=cross_length, lineWidth=8, size=10, closeShape=False, lineColor='white')

# ============================================================
# TEXT MESSAGES
# ============================================================
# Welcome / general instructions
instr = '''
Welcome, and thank you for participating in this experiment!\n
The experiment will be composed of a number of blocks, each containing one out of two tasks:
In one task "HEART TASK", you will be asked to close your eyes, focus on your heart and count your heartbeats for a number of trials.
Each trial will begin and end with a clear tone, and after each trial you will be asked to report how many heartbeats you counted.
In the other task "TONE TASK", you will be asked to close your eyes, focus on the faint tones being played and count the number of tones (not counting the clear start and stop tones).
After each trial, you will be asked to report how many tones/heartbeats you counted (depending on the task), and how confident you are in your report.
Each block will run for 10 trials, and there will be a total of 12 blocks.
(Press any key to begin)
'''

brk = '''Congrats! You've just completed a block. Please take a moment to relax, and press any button to start the next task'''

end = '''The experiment is over, thank you very much for participating!
Please, call the experimenter before leaving.
'''

# Accuracy questions
acc_int = '''How many heartbeats did you count?

(Note: You need to press any number key once before you can type your response)




\n'''
acc_ext = '''How many tones did you count?

(Note: You need to press any number key once before you can type your response)




\n'''

# Alert screens
alert_int = '''HEART TASK

Close your eyes'''
alert_ext = '''TONE TASK

Close your eyes'''

# Task introduction instructions
intro_int = '''HEART TASK

In the following task, you will be asked to close your eyes, and focus your attention on your heart and try to feel your heartbeats.
When you hear the clear tone, start counting your heartbeats. When you hear the clear tone again, stop counting and open your eyes, and report how many heartbeats you felt.
After that, you will be asked to rate how confident you are about your response.
(Press any key to start)
'''
intro_ext = '''TONE TASK

In the following task, you will be asked to close your eyes, and focus your attention on the faint tones being played.
When you hear the clear tone, start counting the faint tones. When you hear the clear tone again, open your eyes and report how many tones you counted (not counting the clear start and stop tones).
After that, you will be asked to rate how confident you are about your response.
(Press any key to start)
'''

# Training instructions
training_intro = '''Before starting the main experiment task, you will now be going through a quick training session.
During the task, you will be asked to close your eyes, listen for a number of faint tones, and count how many are played.
You will hear a loud and clear tone marking the start of the trial. After you hear this, listen for the faint tones and count how many you hear. Please keep your eyes closed.
When you hear the clear tone again, the trial is over. You may then open your eyes, and report how many faint tones you counted (NOT including the clear start and stop tones), and rate how confident you are in your response.
This will be repeated for a few trials.
(Press any button to start)
'''

training_complete = '''Well done! The training is now over.
Now we move on to the main task
'''

# PAS (Perceptual Awareness Scale) - Interoceptive
PAS_int_text = '''How confident are you in your report?\n
1 - I did not sense my heartbeats; I am completely guessing about the number of beats\n
2 - I sensed something about my heart, but I had no idea what I was counting, and I have no confidence at all in my counting\n
3 - I sporadically or faintly picked up on my heart beat; my counting is based on something, but it may be off by a small margin\n
4 - I clearly sensed my heart beat, and have full confidence in my count
'''
PAS_int_text1 = '''
1 - I did not sense my heartbeats; I am completely guessing about the number of beats\n
'''
PAS_int_text2 = '''
2 - I sensed something about my heart, but I had no idea what I was counting, and I have no confidence at all in my counting\n
'''
PAS_int_text3 = '''
3 - I sporadically or faintly picked up on my heart beat; my counting is based on something, but it may be off by a small margin\n
'''
PAS_int_text4 = '''
4 - I clearly sensed my heart beat, and have full confidence in my count
'''

# PAS - Exteroceptive
PAS_ext_text = '''How confident are you in your report?\n
1 - I did not hear any tones; I am completely guessing about the number of tones\n
2 - I heard something like tones, but I had no idea what I was counting, and I have no confidence at all in my counting\n
3 - I sporadically or faintly picked up on the tones; my counting is based on something, but it may be off by a small margin\n
4 - I clearly heard the tones, and have full confidence in my counting
'''
PAS_ext_text1 = '''
1 - I did not hear any tones; I am completely guessing about the number of tones\n
'''
PAS_ext_text2 = '''
2 - I heard something like tones, but I had no idea what I was counting, and I have no confidence at all in my counting\n
'''
PAS_ext_text3 = '''
3 - I sporadically or faintly picked up on the tones; my counting is based on something, but it may be off by a small margin\n
'''
PAS_ext_text4 = '''
4 - I clearly heard the tones, and have full confidence in my counting
'''

# ============================================================
# VISUAL STIMULUS OBJECTS
# ============================================================
startup_text = visual.TextStim(win, text='Program starting. Press any key to begin.', height=msgheight, pos=msgpos)
instr_text = visual.TextStim(win, text=instr, height=msgheight, pos=msgpos)
instr_int_text = visual.TextStim(win, text=intro_int, height=msgheight, pos=msgpos)
instr_ext_text = visual.TextStim(win, text=intro_ext, height=msgheight, pos=msgpos)
break_text = visual.TextStim(win, text=brk, height=msgheight, pos=msgpos)
end_text = visual.TextStim(win, text=end, height=msgheight, pos=msgpos)
acc_int_text = visual.TextStim(win, text=acc_int, height=msgheight, pos=msgpos)
acc_ext_text = visual.TextStim(win, text=acc_ext, height=msgheight, pos=msgpos)
alert_int_text = visual.TextStim(win, text=alert_int, height=(msgheight * 2), color='red', pos=msgpos)
alert_ext_text = visual.TextStim(win, text=alert_ext, height=(msgheight * 2), color='blue', pos=msgpos)
error_warning = visual.TextStim(win, text='Invalid answer. Please try again.', height=msgheight, pos=msgpos)
training_instr_text = visual.TextStim(win, text=training_intro, height=msgheight, pos=msgpos)
training_complete_text = visual.TextStim(win, text=training_complete, height=msgheight, pos=msgpos)
CapturedResponseString = visual.TextStim(win, text='', height=msgheight, pos=(0, -15))

# PAS visual stimuli
PAS_int_conf = visual.TextStim(win, text=PAS_int_text, height=msgheight, pos=msgpos)
PAS_int_conf1 = visual.TextStim(win, text=PAS_int_text1, height=msgheight, pos=(0, 175), color='red')
PAS_int_conf2 = visual.TextStim(win, text=PAS_int_text2, height=msgheight, pos=(0, 18), color='red')
PAS_int_conf3 = visual.TextStim(win, text=PAS_int_text3, height=msgheight, pos=(0, -157), color='red')
PAS_int_conf4 = visual.TextStim(win, text=PAS_int_text4, height=msgheight, pos=(0, -280), color='red')

PAS_ext_conf = visual.TextStim(win, text=PAS_ext_text, height=msgheight, pos=msgpos)
PAS_ext_conf1 = visual.TextStim(win, text=PAS_ext_text1, height=msgheight, pos=(0, 175), color='red')
PAS_ext_conf2 = visual.TextStim(win, text=PAS_ext_text2, height=msgheight, pos=(0, 18), color='red')
PAS_ext_conf3 = visual.TextStim(win, text=PAS_ext_text3, height=msgheight, pos=(0, -157), color='red')
PAS_ext_conf4 = visual.TextStim(win, text=PAS_ext_text4, height=msgheight, pos=(0, -280), color='red')

# ============================================================
# HELPER FUNCTIONS
# ============================================================
# Mapping from PsychoPy key names to numeric characters
NUM_KEY_MAP = {
    'num_0': '0', '0': '0', 'num_1': '1', '1': '1',
    'num_2': '2', '2': '2', 'num_3': '3', '3': '3',
    'num_4': '4', '4': '4', 'num_5': '5', '5': '5',
    'num_6': '6', '6': '6', 'num_7': '7', '7': '7',
    'num_8': '8', '8': '8', 'num_9': '9', '9': '9',
}
SPECIAL_KEY_MAP = {'space': ' ', 'period': '.', 'comma': ','}
IGNORE_KEYS = {'lshift', 'rshift'}

captured_string = ''


def cleanup_and_quit():
    """Close all result files and exit the experiment."""
    results_training.close()
    results_int.close()
    results_ext.close()
    win.close()
    core.quit()


def update_response_display(captured_string):
    """Update the on-screen text as the participant types."""
    CapturedResponseString.setText(captured_string)
    CapturedResponseString.draw()
    acc_text.draw()
    win.flip()


def collect_numeric_input():
    """Collect a numeric string from the participant via keyboard.

    Waits for the participant to type digits and press return.
    Validates that the input contains only digits before accepting.
    Returns the captured string.
    """
    global captured_string
    subject_response_finished = 0

    while subject_response_finished == 0:
        for key in event.getKeys():
            if key in ['escape']:
                cleanup_and_quit()
            elif key in ['return']:
                if captured_string.isdigit():
                    subject_response_finished = 1
                else:
                    error_warning.draw()
                    win.flip()
                    captured_string = ''
            elif key in ['delete', 'backspace']:
                captured_string = captured_string[:-1]
                update_response_display(captured_string)
            elif key in NUM_KEY_MAP:
                captured_string = captured_string + NUM_KEY_MAP[key]
                update_response_display(captured_string)
            elif key in SPECIAL_KEY_MAP:
                captured_string = captured_string + SPECIAL_KEY_MAP[key]
                update_response_display(captured_string)
            elif key in IGNORE_KEYS:
                pass
            else:
                captured_string = captured_string + key
                update_response_display(captured_string)

    return captured_string


def collect_PAS_response():
    """Collect a PAS (Perceptual Awareness Scale) rating (1-4) from the participant.

    Displays the PAS scale and waits for the participant to select 1-4 and press return.
    Returns the selected rating as a string.
    """
    PAS_response = ''
    subject_response_finished = 0

    PAS_conf.draw()
    win.flip()

    while subject_response_finished == 0:
        for key in event.getKeys():
            if key in ['escape']:
                cleanup_and_quit()
            elif key in ['return'] and PAS_response in ['1', '2', '3', '4']:
                subject_response_finished = 1
            elif key in ['num_1', '1']:
                PAS_response = '1'
                PAS_conf.draw()
                PAS_conf1.draw()
                win.flip()
            elif key in ['num_2', '2']:
                PAS_response = '2'
                PAS_conf.draw()
                PAS_conf2.draw()
                win.flip()
            elif key in ['num_3', '3']:
                PAS_response = '3'
                PAS_conf.draw()
                PAS_conf3.draw()
                win.flip()
            elif key in ['num_4', '4']:
                PAS_response = '4'
                PAS_conf.draw()
                PAS_conf4.draw()
                win.flip()

    return PAS_response


def get_next_interval(n_intervals):
    """Get the next trial interval from the intervals list, reshuffling if needed.

    Returns (trialLength, updated n_intervals).
    """
    try:
        trialLength = intervals[n_intervals]
        n_intervals += 1
    except IndexError:
        random.shuffle(intervals)
        n_intervals = 0
        trialLength = intervals[n_intervals]
        n_intervals += 1
    return trialLength, n_intervals


# ============================================================
# TONE STIMULUS PARAMETERS
# ============================================================
# Staircase volume settings
tone_vol_low = 0.005
tone_vol_high = 0.05
tone_vol = tone_vol_high
vol_step = 0.0005

# Start/stop tick sound
tick = sound.Sound(600, secs=0.1, sampleRate=44100, loops=0)

# Tone stimulus parameters
tone_Hz = 250
tone_secs = 0.1
tone_sampleRate = 44100

toneStim = sound.Sound(value=tone_Hz, secs=tone_secs, sampleRate=tone_sampleRate, loops=0)
toneStim.setVolume(tone_vol)

# Inter-stimulus interval (ISI) options
toneStim_isi_list1 = [0.5, 0.75, 1, 1.25, 1.5]
toneStim_isi_list2 = [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5]

# Tone pitch randomization range (Hz)
tone_Hz_range_min = 225
tone_Hz_range_max = 275

toneStim_isi_lists = [toneStim_isi_list1, toneStim_isi_list2]
toneStim_isi_list = random.choice(toneStim_isi_lists)
toneStim_isi = random.choice(toneStim_isi_list)

# Target hit rate range for staircase
hitRange_min = 0.55
hitRange_max = 0.85

# ============================================================
# STARTUP SCREEN
# ============================================================
startup_text.draw()
win.flip()
event.waitKeys()

# ============================================================
# TRAINING PHASE (Staircase volume calibration)
# ============================================================
training_instr_text.draw()
win.flip()
event.waitKeys()
fixation.draw()
win.flip()
core.wait(1)

random.shuffle(intervals)

# Training uses exteroceptive (tone) task only
alert_text = alert_ext_text
acc_text = acc_ext_text
instr_ext_text.draw()
win.flip()
event.waitKeys()

# Set PAS scale to exteroceptive version for training
PAS_conf = PAS_ext_conf
PAS_conf1 = PAS_ext_conf1
PAS_conf2 = PAS_ext_conf2
PAS_conf3 = PAS_ext_conf3
PAS_conf4 = PAS_ext_conf4

PAS_response = ''
n_intervals = 0
n_testTrials = 0

# Consecutive-hits counters for staircase convergence
volHigh_hitCounter = 0
volLow_hitCounter = 0

# Volume randomizer state: alternates between high and low volume
volHighLow_current = 0   # 0 = high volume, 1 = low volume
volHighLow_str = 'High'
volHighLow_pick = 0
volHighLow_consecutive = 1  # Limit to max 2 consecutive trials at same volume

if not training_running:
    vol_final = (tone_vol_high + tone_vol_low) / 2  # Default if training is skipped

while training_running:
    n_testTrials += 1
    print('trial nr ' + str(n_testTrials))

    # Repick ISI list
    toneStim_isi_lists = [toneStim_isi_list1, toneStim_isi_list2]
    toneStim_isi_list = random.choice(toneStim_isi_lists)
    toneStim_isi = random.choice(toneStim_isi_list)

    alert_text.draw()
    win.flip()

    # Set stimulus volume based on current high/low staircase arm
    if volHighLow_current == 0:
        tone_vol = tone_vol_high
        volHighLow_str = 'High'
    elif volHighLow_current == 1:
        tone_vol = tone_vol_low
        volHighLow_str = 'Low'

    # Get trial length from interval list
    trialLength, n_intervals = get_next_interval(n_intervals)

    core.wait(iti)

    tick = sound.Sound(600, secs=0.1, sampleRate=44100, loops=0)
    tick.play()  # START tone
    core.wait(0.1)

    if ppPresent:
        sendParallelTrigger(pport_address_out, triggerStartTraining)

    startTime = time.time()
    stimTime = startTime + 1
    toneStim_check = False
    toneStim_counter = 0

    # Trial loop: generate tones at random intervals
    while time.time() < (startTime + trialLength):
        if time.time() > (startTime + 1):  # 1-second buffer after start
            if time.time() > (stimTime + toneStim_isi) and toneStim_check == False:
                tone_Hz = random.randint(tone_Hz_range_min, tone_Hz_range_max)
                toneStim = sound.Sound(value=tone_Hz, secs=tone_secs, sampleRate=tone_sampleRate, loops=0, blockSize=256, hamming=False)
                toneStim.setVolume(tone_vol)
                toneStim.play()
                if ppPresent:
                    sendParallelTrigger(pport_address_out, mrkrCode_tone)
                toneStim_counter += 1
                stimTime = time.time()
                toneStim_isi = random.choice(toneStim_isi_list)
                toneStim_check = True
            elif toneStim_check == True and time.time() > (stimTime + toneStim_isi):
                toneStim_check = False
            else:
                core.wait(max(0, stimTime + toneStim_isi - time.time()))

    tick = sound.Sound(600, secs=0.1, sampleRate=44100, loops=0)
    tick.play()  # END tone
    core.wait(0.1)

    win.flip()
    if ppPresent:
        sendParallelTrigger(pport_address_out, triggerEndTraining)

    core.wait(0.5)

    # Collect participant's count
    acc_text.draw()
    win.flip()
    event.waitKeys()  # First keypress activates input mode

    captured_string = ''
    collect_numeric_input()

    # Collect PAS rating
    PAS_response = collect_PAS_response()

    # Calculate hit rate
    tones_reported = int(captured_string)
    if toneStim_counter == 0:
        toneStim_counter = 1  # Avoid ZeroDivisionError
    trial_hitRate = 1 - (abs(float(toneStim_counter - tones_reported)) / float(toneStim_counter))

    print('reported: ' + str(tones_reported))
    print('actual: ' + str(toneStim_counter))
    print('hit rate: ' + str(trial_hitRate))
    print('volume: ' + str(tone_vol))
    print('confidence: ' + str(PAS_response))

    # Write training results
    results_training.write(
        str(n_testTrials) + ',' +
        str(tones_reported) + ',' +
        str(toneStim_counter) + ',' +
        str(trial_hitRate) + ',' +
        str(tone_vol) + ',' +
        str(volHighLow_str) + ',' +
        str(PAS_response) + '\n')
    results_training.flush()

    # Adjust volume step size
    if n_testTrials < 10:
        vol_step = tone_vol * 0.5
    else:
        vol_step = tone_vol * 0.25  # Finer adjustments after 10 trials

    # Adjust volume based on hit rate
    if trial_hitRate < hitRange_min:
        print('vol increased')
        tone_vol += vol_step
    elif trial_hitRate > hitRange_max and trial_hitRate <= 1:
        print('vol decreased')
        tone_vol -= vol_step
    elif trial_hitRate > 1:
        print('overshoot')
    else:
        print('vol the same')
    print('######################')

    # Track consecutive hits within target range
    if hitRange_min <= trial_hitRate <= hitRange_max:
        if volHighLow_current == 0:
            volHigh_hitCounter += 1
        elif volHighLow_current == 1:
            volLow_hitCounter += 1
    else:
        if volHighLow_current == 0:
            volHigh_hitCounter = 0
        elif volHighLow_current == 1:
            volLow_hitCounter = 0

    # Save current volume to the appropriate staircase arm
    if volHighLow_pick == 0:
        tone_vol_high = tone_vol
    elif volHighLow_pick == 1:
        tone_vol_low = tone_vol

    toneStim_isi_list = random.choice(toneStim_isi_lists)

    # Randomly pick high or low volume for next trial (max 2 consecutive of same)
    volHighLow_pick = random.randint(0, 1)
    if volHighLow_pick == volHighLow_current:
        volHighLow_consecutive += 1
        if volHighLow_consecutive > 2:
            volHighLow_pick = 1 - volHighLow_pick  # Flip
            volHighLow_consecutive = 1
    volHighLow_current = volHighLow_pick
    if volHighLow_current == 0:
        tone_vol = tone_vol_high
    elif volHighLow_current == 1:
        tone_vol = tone_vol_low

    # Reset for next trial
    captured_string = ''
    PAS_response = ''

    # Check staircase convergence
    if n_testTrials >= 10:
        if volHigh_hitCounter >= 2 or volLow_hitCounter >= 2:
            vol_final = tone_vol
            training_running = False
    if n_testTrials >= 20:
        # Fallback: use midpoint if no convergence after 20 trials
        vol_final = (tone_vol_high + tone_vol_low) / 2
        training_running = False

tone_vol = vol_final

training_complete_text.draw()
win.flip()
core.wait(5)

# ============================================================
# MAIN TASK
# ============================================================
instr_text.draw()
win.flip()
core.wait(3)
event.waitKeys()

# Randomly choose starting condition, then alternate
condition = random.choice(conditions)

count_b = 0        # Block counter
count_t_tot = 0    # Total trial counter
HB_tracker = []    # Track heartbeat counts across interoceptive trials

tone_vol = vol_final  # Use calibrated volume from training

for b in n_blocks:

    tracker_hitRate = []  # Per-block hit rate tracker
    count_b += 1
    count_t_cond = 0
    n_intervals = 0
    hitRate_previous = 0

    # Set up condition-specific stimuli and PAS scale
    if condition == 0:  # Interoceptive (heart) task
        alert_text = alert_int_text
        acc_text = acc_int_text
        instr_int_text.draw()
        PAS_conf = PAS_int_conf
        PAS_conf1 = PAS_int_conf1
        PAS_conf2 = PAS_int_conf2
        PAS_conf3 = PAS_int_conf3
        PAS_conf4 = PAS_int_conf4
    elif condition == 1:  # Exteroceptive (tone) task
        alert_text = alert_ext_text
        acc_text = acc_ext_text
        instr_ext_text.draw()
        PAS_conf = PAS_ext_conf
        PAS_conf1 = PAS_ext_conf1
        PAS_conf2 = PAS_ext_conf2
        PAS_conf3 = PAS_ext_conf3
        PAS_conf4 = PAS_ext_conf4

    win.flip()
    core.wait(2)
    event.waitKeys()

    random.shuffle(intervals)

    # Trial loop
    for i in range(n_trials):

        count_t_tot += 1
        count_t_cond += 1
        alert_text.draw()
        win.flip()

        print('block ' + str(count_b) + ', trial ' + str(count_t_cond) + ', cond ' + str(condition))

        # Get trial length
        trialLength, n_intervals = get_next_interval(n_intervals)

        # Randomly pick ISI range for this trial
        toneStim_isi_list = random.choice(toneStim_isi_lists)
        toneStim_isi = random.choice(toneStim_isi_list)

        core.wait(iti)

        # Wait for heartbeat before starting trial (sync to cardiac cycle)
        start_beat = False
        if ppPresent:
            while start_beat == False:
                signal_in = readParallelTrigger(pport_address_in)
                if signal_in == heartbeatCode:
                    start_beat = True

        if ppPresent:
            sendParallelTrigger(pport_address_out, triggerStart)

        tick = sound.Sound(600, secs=0.1, sampleRate=44100)
        tick.play()  # START tone
        core.wait(0.1)

        startTime = time.time()
        stimTime = startTime + 2
        toneStim_check = False
        toneStim_counter = 0
        HB_counter = 0
        beat_1st = False
        beat_time_list = []

        # Main trial loop
        while time.time() < (startTime + trialLength + 0.1):
            if condition == 0:  # Interoceptive: count heartbeats via parallel port
                if ppPresent:
                    signal_in = readParallelTrigger(pport_address_in)
                    if signal_in == heartbeatCode and beat_1st == True and time.time() > beat_time + 0.1:
                        HB_counter += 1
                        beat_time = time.time()
                        beat_time_list.append(beat_time)
                    elif signal_in == heartbeatCode and beat_1st == False:
                        beat_1st = True
                        beat_time = time.time()
                        beat_time_list.append(beat_time)

            elif condition == 1:  # Exteroceptive: play tones at random intervals
                if time.time() > (startTime + 1) and time.time() < (startTime + trialLength - 0.8):
                    if time.time() > (stimTime + toneStim_isi) and toneStim_check == False:
                        tone_Hz = random.randint(tone_Hz_range_min, tone_Hz_range_max)
                        toneStim = sound.Sound(value=tone_Hz, secs=tone_secs, sampleRate=tone_sampleRate)
                        toneStim.setVolume(tone_vol)
                        toneStim.play()
                        if ppPresent:
                            sendParallelTrigger(pport_address_out, mrkrCode_tone)
                        toneStim_counter += 1
                        stimTime = time.time()
                        toneStim_isi = random.choice(toneStim_isi_list)
                        toneStim_check = True
                    elif toneStim_check == True and time.time() > (stimTime + toneStim_isi):
                        toneStim_check = False
                    else:
                        core.wait(max(0, stimTime + toneStim_isi - time.time()))

        tick = sound.Sound(600, secs=0.1, sampleRate=44100)
        tick.play()  # END tone
        core.wait(0.1)
        win.flip()

        if ppPresent:
            sendParallelTrigger(pport_address_out, triggerEnd)

        core.wait(0.5)

        # Collect participant's count
        acc_text.draw()
        win.flip()
        event.waitKeys()  # First keypress activates input mode

        captured_string = ''
        collect_numeric_input()

        tone_vol_current = tone_vol  # Save before potential staircase adjustment
        print('current volume: ' + str(tone_vol_current))

        # Calculate hit rate (accuracy)
        stim_reported = int(captured_string)
        if HB_counter == 0 and condition == 0:
            if len(HB_tracker) > 0:
                HB_counter = int(round(np.mean(HB_tracker)))
            else:
                HB_counter = 1  # Fallback for first trial
        if toneStim_counter == 0:
            toneStim_counter = 1  # Avoid ZeroDivisionError
        if condition == 0:
            HB_tracker.append(HB_counter)
            trial_hitRate = 1 - (abs(float(HB_counter - stim_reported)) / float(HB_counter))
            print("number of reported heartbeats: " + str(stim_reported))
            print("number of caught heartbeats via RecView: " + str(HB_counter))
        elif condition == 1:
            trial_hitRate = 1 - (abs(float(toneStim_counter - stim_reported)) / float(toneStim_counter))
            print("number of reported tones: " + str(stim_reported))
            print("number of presented tones: " + str(toneStim_counter))
        print('trial hit-rate: ' + str(trial_hitRate))

        # Calculate beats per minute
        bpm = (float(HB_counter) / float(trialLength)) * 60

        # Staircase: track hit-rate and update difficulty
        vol_step = tone_vol * 0.25
        tracker_hitRate.append(trial_hitRate)

        if condition == 0:  # Interoceptive: update hit-rate range at end of block
            if i == max(range(n_trials)):
                avg_hitRate = np.mean(tracker_hitRate)
                hitRange_min = avg_hitRate - 0.15
                hitRange_max = avg_hitRate + 0.15
                if hitRange_min < 0:
                    hitRange_min = 0
                elif hitRange_min > 0.85:
                    hitRange_min = 0.85
                if hitRange_max > 1:
                    hitRange_max = 1
                elif hitRange_max < 0.15:
                    hitRange_max = 0.15
                print('New hit-rate range: ' + str(hitRange_min) + '-' + str(hitRange_max))

        elif condition == 1:  # Exteroceptive: adjust volume if 2 consecutive misses
            if trial_hitRate < hitRange_min:
                hitRate_current = 1  # Low
                print('hit-rate low')
            elif trial_hitRate > hitRange_max and trial_hitRate <= 1:
                hitRate_current = 2  # High
                print('hit-rate high')
            elif trial_hitRate > 1:
                hitRate_current = 0  # Overshoot
                print('hit-rate overshoot')
            else:
                hitRate_current = 0  # Within range
                print('hit-rate within range')

            if hitRate_current == hitRate_previous and hitRate_current != 0:
                if hitRate_current == 1:
                    tone_vol += vol_step
                elif hitRate_current == 2:
                    tone_vol -= vol_step
            hitRate_previous = hitRate_current

        # Collect PAS rating
        PAS_response = collect_PAS_response()
        print('PAS rating: ' + str(PAS_response))

        # Save behavioural results
        if condition == 0:
            hitRate_int = trial_hitRate
            results_int.write(
                str(count_t_tot) + ',' +
                str(count_t_cond) + ',' +
                str(trialLength) + ',' +
                str(count_b) + ',' +
                str(condition) + ',' +
                str(HB_counter) + ',' +
                str(captured_string) + ',' +
                str(hitRate_int) + ',' +
                str(PAS_response) + ',' +
                str(bpm) + '\n')
            results_int.flush()

        elif condition == 1:
            hitRate_ext = trial_hitRate
            results_ext.write(
                str(count_t_tot) + ',' +
                str(count_t_cond) + ',' +
                str(trialLength) + ',' +
                str(count_b) + ',' +
                str(condition) + ',' +
                str(tone_vol_current) + ',' +
                str(toneStim_counter) + ',' +
                str(captured_string) + ',' +
                str(hitRate_ext) + ',' +
                str(PAS_response) + ',' +
                str(bpm) + '\n')
            results_ext.flush()

        # Reset for next trial
        captured_string = ''
        PAS_response = ''
        print('#############')

    # Alternate condition at end of each block
    if condition == 0:
        condition = 1
    elif condition == 1:
        condition = 0

    if count_b == (max(n_blocks) + 1):  # Last block: show end screen
        end_text.draw()
        win.flip()
        core.wait(3)
        break
    else:  # Between blocks: show break screen
        break_text.draw()
        win.flip()
        event.waitKeys()

# ============================================================
# CLEANUP
# ============================================================
results_int.close()
results_ext.close()
results_training.close()

print('data saved')

core.wait(final_t)
win.close()
