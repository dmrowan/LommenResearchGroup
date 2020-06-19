timetab['timedif']= timetab['STOP'] - timetab['START']
totaltime=0
startime = timetab['START']
for counter in range(np.size(timetab[0])):
    totaltime += timetab['timedif'][counter]
    if (totaltime >= timewidth):
        diff = totaltime-timewidth
        totaltime -= diff
        endtime = timetab['STOP'][counter] - diff
        # make profiles
        goodtimes = np.where(time>startime && time < endtime)
        make the profile
        fit it
        write it to a file
        startime=endtime
        totaltime=0
