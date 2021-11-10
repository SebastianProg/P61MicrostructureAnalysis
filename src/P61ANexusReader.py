import h5py
import numpy as np
import pandas as pd
import os

import basics.functions as gf
import filehandling.general as gfh
import filehandling.specific as sfh


# implementation by Gleb
class P61ANexusReader:
    # linear energy calibration for ch0 and ch1 (but given as quadratic coefficients)
    calib = np.array([[0, 0.0505, -0.1], [0, 0.0499999999999995, -0.0249999999990678]])
    # linear dead time calibration for ch0 and ch1 only based on first Fe measurements with different graphite plates
    # dtcalib = np.array([[2.50101945306792E-06, 0.0163784591685192], [0.0000056494547880336, 0.0379388442097675]])
    # quadratic dead time calibration for ch0 and ch1 based on both Fe measurements with different graphite plates
    dtcalib = np.array([[2.6927340476936e-14, 2.48127889225408e-06, 0.0224972569377339],
                        [-2.21978987294787e-13, 5.96715902282922e-06, -0.0320589089676879]])
    ch0 = ('entry', 'instrument', 'xspress3', 'channel00')
    ch1 = ('entry', 'instrument', 'xspress3', 'channel01')
    all_event = ('scaler', 'allevent')
    all_good = ('scaler', 'allgood')
    time = ('scaler', 'time')
    reset_ticks = ('scaler', 'resetticks')
    hist = ('histogram', )
    columns = ('DataX', 'DataY', 'DataID', 'Channel', 'ScreenName', 'Active', 'Color', 'DeadTime', 'CountTime')

    def __init__(self):
        self.q_app = None

    def validate(self, f_name):
        hists = False
        try:
            with h5py.File(f_name, 'r') as f:
                hists |= '/'.join(self.ch0 + self.hist) in f
                hists |= '/'.join(self.ch1 + self.hist) in f
        except Exception:
            return False

        return hists

    def read(self, f_name, sum_frames=True):
        # kev_per_bin = 5E-2
        if self.q_app is not None:
            result = pd.DataFrame(columns=self.q_app.data.columns)
        else:
            result = pd.DataFrame(columns=self.columns)

        try:
            with h5py.File(f_name, 'r') as f:
                for ii, channel in enumerate((self.ch0, self.ch1)):
                    if '/'.join(channel + self.hist) not in f:
                        continue

                    frames = np.sum(f['/'.join(channel + self.hist)], axis=0)
                    frames[:20] = 0.0
                    frames[-1] = 0.0
                    # calculation of energies
                    # kev = np.arange(frames.shape[0]) * kev_per_bin
                    # if ii == 0:
                    #     kev = np.arange(frames.shape[0]) * 0.050494483569344 + 0.029899315869827
                    # elif ii == 1:
                    #     kev = np.arange(frames.shape[0]) * 0.04995786201326 + 0.106286326963684
                    # else:
                    #     kev = (np.arange(frames.shape[0]) + 0.5) * kev_per_bin
                    ch_vals = np.arange(frames.shape[0]) + 1
                    kev = self.calib[ii, 0] * ch_vals ** 2 + self.calib[ii, 1] * ch_vals + self.calib[ii, 2]

                    if self.q_app is not None:
                        row = {c: None for c in self.q_app.data.columns}
                    else:
                        row = {c: None for c in self.columns}

                    # if ('/'.join(channel + self.all_event) in f) and ('/'.join(channel + self.all_good) in f):
                    #     allevent = np.sum(f['/'.join(channel + self.all_event)], axis=0)
                    #     allgood = np.sum(f['/'.join(channel + self.all_good)], axis=0)
                    #     row.update({'DeadTime': (1. - allgood / allevent) * 100})

                    count_time = np.sum(f['/'.join(channel + self.time)], axis=0) * 1.25e-8
                    resetticks = np.sum(f['/'.join(channel + self.reset_ticks)], axis=0)
                    dt_x = resetticks / count_time

                    row.update({
                        'DataX': kev,
                        'DataY': frames,
                        'DataID': f_name + ':' + '/'.join(channel),
                        'Channel': ii,
                        'ScreenName': os.path.basename(f_name) + ':' + '%02d' % ii,
                        'Active': True,
                        'CountTime': count_time,
                        'DeadTime': self.dtcalib[ii, 0] * dt_x ** 2 + self.dtcalib[ii, 1] * dt_x + self.dtcalib[ii, 2]
                    })

                    if self.q_app is not None:
                        row.update({'Color': next(self.q_app.params['ColorWheel'])})
                    result.loc[result.shape[0]] = row
        except Exception:
            print('Problem with file: ' + f_name)
            return False

        result = result.astype('object')
        result[pd.isna(result)] = None
        return result

    def main(self):
        differentOutputFolder = True
        # differentOutputFolder = False
        # select nexus files
        fileNames = gfh.requestFiles((("Nexus files", "*.nxs"),), "Select P61A file", "on")
        pathName, _ = gfh.fileparts(fileNames[0])
        # read all data from selected files
        reader = P61ANexusReader()
        data = pd.DataFrame(columns=reader.columns)
        for f_name in fileNames:
            data = pd.concat((data, reader.read(f_name)), ignore_index=True)
        # select output folder if wanted
        if differentOutputFolder:
            outputPath = gfh.requestDirectory(dialogTitle='Select output folder', folder=pathName)
        else:
            outputPath = pathName
        # export data to spectrum files
        for i in range(data.shape[0]):
            curData = data.iloc[i]
            header = ['SCREEN_NAME=' + curData['ScreenName'], 'CHANNEL=' + str(curData['Channel']),
                'DEADTIME=' + str(curData['DeadTime']), 'TREAL=' + str(curData['CountTime']),
                'SPECTRTXT=' + str(len(curData['DataX']))]
            vals = np.transpose([curData['DataX'], curData['DataY']])
            sfh.writeSpectrumFile(outputPath + '/' + gf.replace(curData['ScreenName'], '.nxs:', '_') + '.txt', vals,
                header)


if __name__ == '__main__':
    P61ANexusReader().main()
