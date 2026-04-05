#!/usr/bin/env python3
"""
external_scorer.py - GenerativeBGCs Validation & Scoring Wrapper

Seamlessly integrates the zero-dependency GenerativeBGCs pipeline with external 
"Gold Standard" third-party algorithms for structural and product validation.

Modes:
  --mode antismash  : Quietly submits .gbk files to antiSMASH API via curl and polls.
  --mode docker     : Uses local Docker to analyze .gbk via DeepBGC machine learning.
  --background      : Daemonizes the polling process so you don't have to wait.
"""

import os
import sys
import glob
import json
import time
import argparse
import subprocess

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
GBK_DIR = os.path.join(RESULTS_DIR, "gbk")
TRACKER_JSON = os.path.join(RESULTS_DIR, "antismash_jobs.json")

class BGCScorer:
    def __init__(self):
        self.jobs = self.load_tracker()

    def load_tracker(self):
        if os.path.exists(TRACKER_JSON):
            try:
                with open(TRACKER_JSON) as f:
                    return json.load(f)
            except:
                pass
        return {}

    def save_tracker(self):
        with open(TRACKER_JSON, 'w') as f:
            json.dump(self.jobs, f, indent=2)

    def run_deepbgc(self, gbk_path):
        """Uses local docker DeepBGC ML model immediately on the file"""
        print(f"[DeepBGC Docker] Instantiating classification on: {os.path.basename(gbk_path)}")
        filename = os.path.basename(gbk_path)
        out_dir = os.path.join(RESULTS_DIR, "deepbgc_scores", filename.replace('.gbk', ''))
        os.makedirs(out_dir, exist_ok=True)

        cmd = [
            "docker", "run", "--rm",
            "-v", f"{RESULTS_DIR}:/data",
            "antibioti/deepbgc", "pipeline",
            f"/data/gbk/{filename}",
            "-o", f"/data/deepbgc_scores/{filename.replace('.gbk', '')}"
        ]
        
        print("  Running: " + " ".join(cmd))
        try:
            # We use Popen so if user Ctrl-C it doesn't break everything
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            for line in proc.stdout:
                line = line.strip()
                if "Score" in line or "Detected" in line or "BGC" in line:
                    print("  [ML Core] " + line)
            proc.wait()
            print(f"  [SUCCESS] Evaluation saved to {out_dir}")
        except FileNotFoundError:
            print("[ERROR] Docker is not installed or not available in PATH.")
        except Exception as e:
            print(f"[ERROR] Docker execution failed: {e}")

    def submit_to_antismash(self, gbk_path, contact="developer@project.local"):
        """Uses raw curl via subprocess to maintain zero-dependency rule"""
        if gbk_path in self.jobs and self.jobs[gbk_path].get("status") in ["running", "completed"]:
            print(f"  [SKIP] {os.path.basename(gbk_path)} already queued: {self.jobs[gbk_path]['job_id']}")
            return

        print(f"[antiSMASH API] Uploading {os.path.basename(gbk_path)} to Leiden server...")
        
        # We rely on system curl to safely send Multipart form data without 'requests' library
        cmd = [
            "curl", "-s", "-X", "POST",
            "https://antismash.secondarymetabolites.org/api/v2.0/upload/",
            "-F", f"data=@{gbk_path}",
            "-F", f"contact={contact}"
        ]
        
        try:
            out = subprocess.check_output(cmd, text=True)
            resp = json.loads(out)
            if "job_id" in resp:
                self.jobs[gbk_path] = {
                    "job_id": resp["job_id"],
                    "status": "queued",
                    "submitted_at": time.time()
                }
                self.save_tracker()
                print(f"  [ACCEPTED] Job ID: {resp['job_id']}")
            else:
                print(f"  [ERROR] Upload failed: {resp}")
        except Exception as e:
            print(f"  [ERROR] Communication failed: {e}")

    def poll_antismash(self):
        """Polls jobs and downloads results."""
        pending = [p for p, data in self.jobs.items() if data.get("status") != "completed"]
        if not pending:
            print("[antiSMASH API] All tracked tasks are already completed.")
            return

        print(f"\n[antiSMASH Tracker] Monitoring {len(pending)} active jobs on server...")
        for p in pending:
            job_id = self.jobs[p]["job_id"]
            url = f"https://antismash.secondarymetabolites.org/api/v2.0/results/{job_id}/"
            try:
                out = subprocess.check_output(["curl", "-s", "-X", "GET", url], text=True)
                resp = json.loads(out)
                
                status = resp.get("status", "unknown")
                print(f"  [Job {job_id}] Remote status: {status}")
                if status == "successful":
                    self.jobs[p]["status"] = "completed"
                    # For now just update status. Downloading full ZIP results is an optional advanced feature
                    print(f"    -> Structure successfully solved! View it at: {url.replace('api/v2.0/results/','upload/')}")
                
            except Exception as e:
                print(f"  [Job {job_id}] Poll failed: {e}")
        self.save_tracker()

def main():
    parser = argparse.ArgumentParser(description="External BGC Scorer (antiSMASH / DeepBGC)")
    parser.add_argument("--mode", choices=["antismash", "docker"], required=True, help="Scoring engine to use.")
    parser.add_argument("--background", action="store_true", help="Daemonize the antiSMASH polling process.")
    args = parser.parse_args()

    scorer = BGCScorer()
    gbk_files = glob.glob(os.path.join(GBK_DIR, "*.gbk"))

    if not gbk_files:
        print("[ERROR] No target .gbk files found in results/gbk")
        sys.exit(1)

    print(f"Found {len(gbk_files)} chimeric templates for Validation.")

    if args.mode == "docker":
        for gbk in gbk_files:
            scorer.run_deepbgc(gbk)
            
    elif args.mode == "antismash":
        for gbk in gbk_files:
            scorer.submit_to_antismash(gbk)
            
        print("\n--- Upload Complete ---")
        if args.background:
            print("[DAEMON] Shifting into background. Will poll server every 5 minutes.")
            pid = os.fork()
            if pid > 0:
                print(f"Daemon PID: {pid}. You can safely close this terminal.")
                sys.exit(0)
            
            # Daemon loop
            with open("/dev/null", "w") as f_null:
                os.dup2(f_null.fileno(), sys.stdout.fileno())
                os.dup2(f_null.fileno(), sys.stderr.fileno())
            
            while True:
                time.sleep(300)
                scorer.poll_antismash()
                if not any(d.get("status") != "completed" for d in scorer.jobs.values()):
                    break
        else:
            scorer.poll_antismash()

if __name__ == "__main__":
    main()
