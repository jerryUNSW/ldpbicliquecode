#!/usr/bin/env python3
"""
Monitor NIPS experiment progress and report status every 2 hours.
"""

import time
import os
import subprocess
import json
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime

def check_file_progress(filename):
    """Check the progress of the NIPS experiment file."""
    if not os.path.exists(filename):
        return "File not found"
    
    try:
        with open(filename, 'r') as f:
            content = f.read()
        
        # Check for completed algorithms
        results = {}
        
        # Check ADV (Algorithm 2)
        if "=== Testing ADV (Algorithm 2) ===" in content:
            if "No-sampling batch results for Q=" in content:
                results['ADV'] = "COMPLETED"
            else:
                results['ADV'] = "RUNNING"
        else:
            results['ADV'] = "NOT STARTED"
        
        # Check ADV+ (Algorithm 3)
        if "=== Testing ADV+ (Algorithm 3) ===" in content:
            # Look for ADV+ specific results line
            lines = content.split('\n')
            adv_plus_started = False
            adv_plus_completed = False
            
            for line in lines:
                if "=== Testing ADV+ (Algorithm 3) ===" in line:
                    adv_plus_started = True
                elif adv_plus_started and "No-sampling batch results for Q=" in line:
                    adv_plus_completed = True
                    break
                elif adv_plus_started and "=== Testing ADV++ (Algorithm 4) ===" in line:
                    break
            
            if adv_plus_completed:
                results['ADV+'] = "COMPLETED"
            elif adv_plus_started:
                results['ADV+'] = "RUNNING"
            else:
                results['ADV+'] = "NOT STARTED"
        else:
            results['ADV+'] = "NOT STARTED"
        
        # Check ADV++ (Algorithm 4)
        if "=== Testing ADV++ (Algorithm 4) ===" in content:
            # Look for ADV++ specific results line
            lines = content.split('\n')
            adv_plus_plus_started = False
            adv_plus_plus_completed = False
            
            for line in lines:
                if "=== Testing ADV++ (Algorithm 4) ===" in line:
                    adv_plus_plus_started = True
                elif adv_plus_plus_started and "No-sampling batch results for Q=" in line:
                    adv_plus_plus_completed = True
                    break
            
            if adv_plus_plus_completed:
                results['ADV++'] = "COMPLETED"
            elif adv_plus_plus_started:
                results['ADV++'] = "RUNNING"
            else:
                results['ADV++'] = "NOT STARTED"
        else:
            results['ADV++'] = "NOT STARTED"
        
        return results
    
    except Exception as e:
        return f"Error reading file: {e}"

def send_email_report(status):
    """Send email report about the progress."""
    try:
        # Load email configuration
        with open('/data1/yizhangh/email_config.json', 'r') as f:
            email_config = json.load(f)
        
        # Get hostname
        hostname = subprocess.check_output(['hostname'], text=True).strip()
        
        # Create email content
        subject = f"NIPS Experiment Progress Report - {hostname}"
        body = f"""
NIPS Experiment Progress Report
Host: {hostname}
Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

Algorithm Status:
- ADV (Algorithm 2): {status.get('ADV', 'UNKNOWN')}
- ADV+ (Algorithm 3): {status.get('ADV+', 'UNKNOWN')}
- ADV++ (Algorithm 4): {status.get('ADV++', 'UNKNOWN')}

File: /data/yizhangh/ldp-pq/testing-nips-no-sampling.txt
        """
        
        # Create email message
        msg = MIMEMultipart()
        msg['From'] = email_config['from_email']
        msg['To'] = 'jerrystat2017@gmail.com'
        msg['Subject'] = subject
        
        msg.attach(MIMEText(body, 'plain'))
        
        # Send email using SMTP
        server = smtplib.SMTP(email_config['smtp_server'], email_config['smtp_port'])
        server.starttls()
        server.login(email_config['from_email'], email_config['password'])
        text = msg.as_string()
        server.sendmail(email_config['from_email'], 'jerrystat2017@gmail.com', text)
        server.quit()
        
        print(f"Email report sent to jerrystat2017@gmail.com at {datetime.now()}")
        
    except Exception as e:
        print(f"Failed to send email: {e}")

def main():
    """Main monitoring loop."""
    filename = "/data/yizhangh/ldp-pq/testing-nips-no-sampling.txt"
    
    print(f"Starting NIPS experiment monitoring...")
    print(f"Monitoring file: {filename}")
    print(f"Reports will be sent every 2 hours")
    print(f"Started at: {datetime.now()}")
    
    while True:
        try:
            # Check current status
            status = check_file_progress(filename)
            
            # Print current status
            print(f"\n=== Status Check at {datetime.now()} ===")
            if isinstance(status, dict):
                for algo, state in status.items():
                    print(f"{algo}: {state}")
            else:
                print(f"Error: {status}")
            
            # Send email report
            if isinstance(status, dict):
                send_email_report(status)
            
            # Wait 2 hours (7200 seconds)
            print(f"Waiting 2 hours until next report...")
            time.sleep(7200)  # 2 hours
            
        except KeyboardInterrupt:
            print("\nMonitoring stopped by user")
            break
        except Exception as e:
            print(f"Error in monitoring loop: {e}")
            time.sleep(300)  # Wait 5 minutes before retrying

if __name__ == "__main__":
    main()
