#!/usr/bin/env python3
"""
Script to send the educational plotting script via email.
"""

import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders
import json
import os

def send_email_with_attachment():
    """Send the educational plotting script via email."""
    try:
        # Load email configuration
        with open('/data1/yizhangh/email_config.json', 'r') as f:
            email_config = json.load(f)
        
        # Email details
        from_email = email_config['from_email']
        to_email = '51255902143@stu.ecnu.edu.cn'
        subject = 'Educational Plotting Script - Publication Quality Plots'
        
        # Create email body
        body = """
Dear Student,

I'm sending you an enhanced version of our comprehensive plotting script with detailed educational comments.

This script demonstrates:
- Professional matplotlib styling for publication-quality figures
- Reading and parsing CSV data files
- Creating line plots with multiple algorithms
- Log scale handling for wide data ranges
- Dynamic subplot layout for multiple datasets
- Exporting high-quality PDF figures

The script includes extensive comments explaining:
- Why each parameter is set (e.g., DPI, font sizes)
- How to customize plots for different use cases
- Best practices for academic figure creation
- Step-by-step workflow explanation

You can use this as a template for creating your own publication-quality plots.

Best regards,
"""
        
        # Create email message
        msg = MIMEMultipart()
        msg['From'] = from_email
        msg['To'] = to_email
        msg['Subject'] = subject
        msg.attach(MIMEText(body, 'plain'))
        
        # Attach the educational plotting script
        filename = 'plot_budget_allocation_styled_educational.py'
        filepath = '/data1/yizhangh/ldp-pq/' + filename
        
        if os.path.exists(filepath):
            with open(filepath, 'rb') as attachment:
                part = MIMEBase('application', 'octet-stream')
                part.set_payload(attachment.read())
            
            encoders.encode_base64(part)
            part.add_header(
                'Content-Disposition',
                f'attachment; filename= {filename}',
            )
            msg.attach(part)
        else:
            print(f"Warning: File {filepath} not found")
            return False
        
        # Send email using SMTP
        server = smtplib.SMTP(email_config['smtp_server'], email_config['smtp_port'])
        server.starttls()
        server.login(from_email, email_config['password'])
        text = msg.as_string()
        server.sendmail(from_email, to_email, text)
        server.quit()
        
        print(f"Email sent successfully to {to_email}")
        print(f"Attachment: {filename}")
        return True
        
    except Exception as e:
        print(f"Failed to send email: {e}")
        return False

if __name__ == '__main__':
    send_email_with_attachment()
