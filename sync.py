# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# import glob
import h5py
import numpy as np


# import os


def parseGNGEvents(events):
    s1s = 30000
    trials = []
    lastTS = -300000
    lastCue = -1
    cueTS = -1
    rsps = -1
    cue = -1

    for eidx in range(len(events)):
        cue = events[eidx][1] & 0x0C
        cueTS = events[eidx][0]
        if cue > 0 and cueTS > (lastTS + s1s):
            if lastCue > 0:
                trials.append([lastTS, lastCue, rsps])
                rsps = -1
            lastCue = cue
            lastTS = cueTS

        if (
                lastCue > 0
                and rsps < 0
                and events[eidx][0] >= lastTS + 30000
                and events[eidx][0] < lastTS + 90000
                and (events[eidx][1] & 0x03) > 0
        ):
            rsps = events[eidx][1] & 0x03

    return trials


def parseDPAEvents(events):
    """
    if no cue,  assume sample
    if has cue history check time diff
    time diff < 3 discard
    time diff > 3 < 10, this cue is test, last cue is sample
    time diff >10

    if sample and test set
        if responsed
            add responsed trial
        else
            add unresponsed trial

        reinit sample test response
        this cue is sample
"""
    s1s = 30000
    trials = []
    lastTS = -300000 * 20
    lastCue = -1
    cueTS = -1
    rsps = -1
    cue = -1
    sample = -1
    test = -1
    sampleTS = -1
    testTS = -1
    Laser = -1
    LaserTS = -1

    for eidx in range(len(events)):
        cue = events[eidx][1] & 0x1E
        cueTS = events[eidx][0]
        
        if cue == 2:
            Laser = 2
            LaserTS = cueTS
        if cue == 16:
            Laser = 1
            LaserTS = cueTS            
        if cue > 2 and cue <16 and cueTS > lastTS + s1s and cueTS < lastTS + 2 * s1s:
            print("error processing evt idx ", eidx)
        elif cue > 2 and cue <16 and cueTS > lastTS + s1s * 2 and cueTS < lastTS + s1s * 8:
            sample = lastCue
            sampleTS = lastTS
            test = cue
            testTS = cueTS

            lastCue = cue
            lastTS = cueTS

        elif cue > 2 and cue <16 and cueTS > lastTS + s1s * 8:
            if sample > 0 and test > 0:
                trials.append(
                    [
                        sampleTS,
                        testTS,
                        np.round(sampleTS / s1s, decimals=3),
                        np.round(testTS / s1s, decimals=3),
                        sample,
                        test,
                        rsps,
                        np.round((testTS - sampleTS) / s1s) - 1,
                        Laser,
                        LaserTS,
                    ]
                )
                sample = -1
                test = -1
                sampleTS = -1
                testTS = -1
                rsps = -1
                Laser = -1
                LaserTS = -1
            lastCue = cue
            lastTS = cueTS

        if (test > 0 and rsps < 0 and events[eidx][0] >= testTS + s1s and events[eidx][0] < lastTS + 2 * s1s and (
                events[eidx][1] & 0x01) > 0):
            rsps = 1
    if sample > 0 and test > 0:
        trials.append(
            [
                sampleTS,
                testTS,
                np.round(sampleTS / s1s, decimals=3),
                np.round(testTS / s1s, decimals=3),
                sample,
                test,
                rsps,
                np.round((testTS - sampleTS) / s1s) - 1,
                Laser,
                LaserTS,
            ]
        )

    return trials

def parseZHADREvents(events):

    """
    if val range from 4 to 24 : sample
    if val == 28 : test
    """

    s1s = 30000
    trials = []
    lastTS = -300000 * 20
    lastCue = -1
    cueTS = -1
    rsps = -1
    cue = -1
    sample = -1
    test = -1
    sampleTS = -1
    testTS = -1

    for eidx in range(len(events)):
        cue = events[eidx][1] & 0x1E
        cueTS = events[eidx][0]
                    
        #if cue == 28 and cueTS > lastTS + 4 * s1s and cueTS < lastTS + 7 * s1s:
        #    test = cue
        #    testTS = cueTS
        #    lastCue = cue
        #    lastTS = cueTS
        
        if cue >= 4 and cue <=25 and cueTS > lastTS + s1s and cueTS < lastTS + 4 * s1s:
            print("error processing evt idx ", eidx)
        elif cue == 28 and cueTS > lastTS + 4 * s1s and cueTS < lastTS + 7 * s1s:
            test = cue
            testTS = cueTS
            sample = lastCue
            sampleTS = lastTS
            lastCue = cue
            lastTS = cueTS
        elif cue >= 4 and cue <= 25 and cueTS > lastTS + 9 * s1s:

            if sample > 0 and test > 0:
                print('here')
                trials.append(
                    [
                        sampleTS,
                        testTS,
                        np.round(sampleTS / s1s, decimals=3),
                        np.round(testTS / s1s, decimals=3),
                        sample,
                        test,
                        rsps,
                        np.round((testTS - sampleTS) / s1s) - 1,
                    ]
                )
                sample = -1
                test = -1
                sampleTS = -1
                testTS = -1
                rsps = -1
            lastCue = cue
            lastTS = cueTS

        if (test > 0 and rsps < 0 and events[eidx][0] >= testTS + s1s and events[eidx][0] < lastTS + 2 * s1s and (
                events[eidx][1] & 0x01) > 0):
            rsps = 1
    if sample > 0 and test > 0:
        trials.append(
            [
                sampleTS,
                testTS,
                np.round(sampleTS / s1s, decimals=3),
                np.round(testTS / s1s, decimals=3),
                sample,
                test,
                rsps,
                np.round((testTS - sampleTS) / s1s) - 1,
            ]
        )

    return trials    




def parseDualTaskEvents(events):

    """
    if no outer cue, assume DPA sample
    if has DPA cue history check time diff
    time diff < 8 discard
    time diff > 8 < 10 this cue is DPA test, last cue is DPA sample
    time diff > 18 abort trial, this cue is DPA sample for next trial, last cue is DPA sample this trial
    1.5s delay response delay
    """
    s1s = 30000
    trials = []
    lastTS = -300000 * 20
    lastCue = -1
    cueTS = -1
    rspsIn = -1
    rspsOut = -1
    cue = -1
    DPASample = -1
    DPATest = -1
    DPASampleTS = -1
    DPATestTS = -1
    innerCue = -1
    lastInnerCue = -1
    innerTS = -1
    lastInnerTS = -300000 * 20
    innerSample = -1
    innerSampleTS = -1
    innerRspCue = -1
    innerRspCueTS = -1

    for eidx in range(len(events)):
        cue = events[eidx][1] & 0x1E
        cueTS = events[eidx][0]
        if cue == 2 or cue == 16:
            innerCue = cue
            innerTS = cueTS
            #print(np.round(innerTS - lastInnerTS) / s1s)
            if innerTS > lastInnerTS + s1s and innerTS < lastInnerTS + 3.5 * s1s:
                innerSample = lastInnerCue
                innerSampleTS = lastInnerTS
                innerRspCue = innerCue
                innerRspCueTS = innerTS

                lastInnerCue = innerCue
                lastInnerTS = innerTS 
            elif innerTS > lastInnerTS + 3.5 * s1s:
                lastInnerCue = innerCue
                lastInnerTS = innerTS
        if (innerRspCue > 0 and rspsIn < 0 and events[eidx][0] >= innerRspCueTS + 0.5 * s1s and events[eidx][0] < lastInnerTS + s1s and (
                events[eidx][1] & 0x01) > 0):
            rspsIn = 1

        if cue > 2 and cue < 16 and cueTS > lastTS + s1s and cueTS < lastTS + 7 * s1s:
            print("error processing evt idx ", eidx)
        elif cue > 2 and cue < 16 and cueTS > lastTS + s1s * 8 and cueTS < lastTS + s1s * 10:
            DPASample = lastCue
            DPASampleTS = lastTS
            DPATest = cue
            DPATestTS = cueTS

            lastCue = cue
            lastTS = cueTS
            print(np.round((DPATestTS - DPASampleTS) / s1s, decimals=3) - 1)
            #and cueTS < lastTS + 15 * s1s
        elif cue > 2 and cue < 16 and cueTS > lastTS + 10 * s1s:
            if DPASample > 0 and DPATest > 0:
                trials.append(
                    [
                        DPASampleTS,
                        DPATestTS,
                        np.round(DPASampleTS / s1s, decimals=3),
                        np.round(DPATestTS / s1s, decimals=3),
                        DPASample,
                        DPATest,
                        rspsOut,
                        np.round((DPATestTS - DPASampleTS) / s1s) - 1,
                        innerSampleTS,
                        innerRspCueTS,
                        np.round(innerSampleTS / s1s, decimals=3),
                        np.round(innerRspCueTS / s1s, decimals=3),
                        innerSample,
                        innerRspCue,
                        rspsIn,
                        np.round((innerRspCueTS - innerSampleTS) / s1s) - 1,
                    ]
                )
                DPASample = -1
                DPATest = -1
                DPASampleTS = -1
                DPATestTS = -1
                rspsOut = -1
                innerSample = -1
                innerRspCue = -1
                innerSampleTS = -1
                innerRspCueTS = -1
                rspsIn = -1
            lastCue = cue
            lastTS = cueTS
            lastInnerCue = innerCue
            lastInnerTS = innerTS
        if (DPATest > 0 and rspsOut < 0 and events[eidx][0] >= DPATestTS + s1s and events[eidx][0] < lastTS + 2 * s1s and (
                events[eidx][1] & 0x01) > 0):
            rspsOut = 1
    if DPASample > 0 and DPATest > 0:
        trials.append(
            [
                DPASampleTS,
                DPATestTS,
                np.round(DPASampleTS / s1s, decimals=3),
                np.round(DPATestTS / s1s, decimals=3),
                DPASample,
                DPATest,
                rspsOut,
                np.round((DPATestTS - DPASampleTS) / s1s) - 1,
                innerSampleTS,
                innerRspCueTS,
                np.round(innerSampleTS / s1s, decimals=3),
                np.round(innerRspCueTS / s1s, decimals=3),
                innerSample,
                innerRspCue,
                rspsIn,
                np.round((innerRspCueTS - innerSampleTS) / s1s) - 1,
            ]
        ) 
    return trials

def getEvents():
    syncs = None
    with h5py.File("sync.hdf5", "r") as fs:
        dset = fs["sync"]
        syncs = np.array(dset, dtype="int8")
        syncs = syncs[0]

    blockCount = 0
    ts = 0
    events = []
    pct = 0
    while ts < (len(syncs))-27:
        if syncs[ts] == 0:
            ts += 1
        else:
            if np.sum(syncs[ts + 6: ts +27])>63:
                state = [ts, 0, 0, 0, 0, 0]
                state[1] = 1 if np.sum(syncs[ts + 7: ts + 8]) > 63 else 0
                state[2] = 1 if np.sum(syncs[ts + 9: ts + 12]) > 127 else 0
                state[3] = 1 if np.sum(syncs[ts + 13: ts + 16]) > 127 else 0
                state[4] = 1 if np.sum(syncs[ts + 17: ts + 20]) > 127 else 0
                state[5] = 1 if np.sum(syncs[ts + 21: ts + 25]) > 127 else 0 
                if (not events) or not np.array_equal(events[-1][1:], state[1:]) :
                    events.append(state)
                ts += 28
            else:
                state = [ts, 0, 0, 0, 0, 0]
                if (not events) or not np.array_equal(events[-1][1:], state[1:]) :
                    events.append(state)
                ts += 6
        # if currPct > pct:
        #     print(currPct)
        #     pct = currPct

        # if syncs[ts] == 0:
        #     blockCount += 1
        #     ts += 1
        # else:
        #     if blockCount < 50:
        #         ts += 1
        #     else:
        #         blockCount = 0
        #         state = [ts, 0, 0, 0, 0, 0]
        #         state[1] = 1 if np.sum(syncs[ts + 5: ts + 7]) > 127 else 0
        #         state[2] = 1 if np.sum(syncs[ts + 9: ts + 14]) > 127 else 0
        #         state[3] = 1 if np.sum(syncs[ts + 14: ts + 16]) > 127 else 0
        #         state[4] = 1 if np.sum(syncs[ts + 16: ts + 20]) > 127 else 0
        #         state[5] = 1 if np.sum(syncs[ts + 20: ts + 25]) > 127 else 0           

        #         # state[6] = 1 if np.sum(syncs[ts + 30: ts + 32]) > 128 else 0
        #         if (not events) or not np.array_equal(events[-1][1:], state[1:]):
        #             events.append(state)
        #         ts += 28
        #         # state = [ts, 0, 0, 0, 0, 0]
        #         # state[1] = 1 if np.sum(syncs[ts + 7: ts + 9]) > 64 else 0
        #         # state[2] = 1 if np.sum(syncs[ts + 10: ts + 11]) > 64 else 0
        #         # state[3] = 1 if np.sum(syncs[ts + 13: ts + 14]) > 64 else 0
        #         # state[4] = 1 if np.sum(syncs[ts + 16: ts + 17]) > 64 else 0
        #         # state[5] = 1 if np.sum(syncs[ts + 20: ts + 21]) > 64 else 0           

        #         # # state[6] = 1 if np.sum(syncs[ts + 30: ts + 32]) > 128 else 0
        #         # if (not events) or not np.array_equal(events[-1][1:], state[1:]):
        #         #     events.append(state)
        #         # if state[5]==0:
        #         #     ts += 20
        #         # else:
        #         #     ts += 24
 
    events = np.array(events)
    return events





def writeEvents(events, trials):
    with h5py.File("events.hdf5", "w") as fw:
        evtDset = fw.create_dataset("events", data=np.array(events, dtype="i4"))
        tDset = fw.create_dataset("trials", data=np.array(trials, dtype="i4"))


def filter_events(events):
    lick_interval = 30000 * 0.05  # 50ms
    cue_interval = 30000 * 0.5
    prev_idx = np.argmax(events[:, 1])
    prev_TS = events[prev_idx, 0]
    for i in range(prev_idx + 1, events.shape[0]):
        if events[i, 1] == 1:
            curr_TS = events[i, 0]
            if curr_TS - prev_TS < lick_interval:
                events[prev_idx:i, 1] = 1
            prev_idx = i
            prev_TS = curr_TS


    prev_idx = np.argmax(events[:, 2])
    prev_TS = events[prev_idx, 0]
    for i in range(prev_idx + 1, events.shape[0]):
        if events[i, 2] == 1:
            curr_TS = events[i, 0]
            if curr_TS - prev_TS < cue_interval:
                events[prev_idx:i, 2] = 1
            prev_idx = i
            prev_TS = curr_TS



    prev_idx = np.argmax(events[:, 3])
    prev_TS = events[prev_idx, 0]
    for i in range(prev_idx + 1, events.shape[0]):
        if events[i, 3] == 1:
            curr_TS = events[i, 0]
            if curr_TS - prev_TS < cue_interval:
                events[prev_idx:i, 3] = 1
            prev_idx = i
            prev_TS = curr_TS

    prev_idx = np.argmax(events[:, 4])
    prev_TS = events[prev_idx, 0]
    for i in range(prev_idx + 1, events.shape[0]):
        if events[i, 4] == 1:
            curr_TS = events[i, 0]
            if curr_TS - prev_TS < cue_interval:
                events[prev_idx:i, 4] = 1
            prev_idx = i
            prev_TS = curr_TS
            
            
    prev_idx = np.argmax(events[:, 5])
    prev_TS = events[prev_idx, 0]
    for i in range(prev_idx + 1, events.shape[0]):
        if events[i, 5] == 1:
            curr_TS = events[i, 0]
            if curr_TS - prev_TS < cue_interval:
                events[prev_idx:i, 5] = 1
            prev_idx = i
            prev_TS = curr_TS
            
    # prev_idx = np.argmax(events[:, 6])
    # prev_TS = events[prev_idx, 0]
    # for i in range(prev_idx + 1, events.shape[0]):
    #     if events[i, 6] == 1:
    #         curr_TS = events[i, 0]
    #         if curr_TS - prev_TS < cue_interval:
    #             events[prev_idx:i, 6] = 1
    #         prev_idx = i
    #         prev_TS = curr_TS
            
    #events[:, 2] = 0
    output = []
    output.append([events[0, 0], events[0, 1] +  events[0, 2] *2 + events[0, 3] * 4 + events[0, 4] * 8+ events[0, 5] * 16])
    for i in range(1, events.shape[0]):
        val = events[i, 1] + events[i, 2] *2 + events[i, 3] * 4 + events[i, 4] * 8 + events[i, 5] * 16 
        if not val == output[-1][1]:
            output.append([events[i, 0], val])
    return output


def runsync():
    events = getEvents()
    events = filter_events(events)
    trials = parseDualTaskEvents(events)
    writeEvents(events, trials)
    return trials


if __name__ == "__main__":
    events = getEvents()
    trials = parseDPAEvents(events)
    writeEvents(events, trials)

# with open('events.csv','w',newline='\r\n') as
