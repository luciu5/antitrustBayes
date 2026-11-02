import subprocess,time,sys,pathlib,csv
command=sys.argv[2:]; output=sys.argv[1]
p=subprocess.Popen(command)
peak_rss=peak_pss=0; scans=0
while p.poll() is None:
    queue=[p.pid]; seen=set(); rss=pss=0
    while queue:
        pid=queue.pop()
        if pid in seen: continue
        seen.add(pid)
        try:
            queue.extend(int(x) for x in pathlib.Path(f'/proc/{pid}/task/{pid}/children').read_text().split())
            for line in pathlib.Path(f'/proc/{pid}/smaps_rollup').read_text().splitlines():
                parts=line.split()
                if parts[0]=='Rss:': rss+=int(parts[1])*1024
                elif parts[0]=='Pss:': pss+=int(parts[1])*1024
        except (FileNotFoundError,ProcessLookupError,PermissionError): pass
    peak_rss=max(peak_rss,rss); peak_pss=max(peak_pss,pss); scans+=1
    time.sleep(.25)
with open(output,'w',newline='') as f:
    w=csv.writer(f);w.writerow(['process_tree_peak_rss_bytes','process_tree_peak_pss_bytes','poll_seconds','scans','exit_status'])
    w.writerow([peak_rss,peak_pss,.25,scans,p.returncode])
sys.exit(p.returncode)
