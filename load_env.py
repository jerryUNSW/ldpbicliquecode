#!/usr/bin/env python3
"""
Load environment variables from project.env file
"""

import os
from pathlib import Path

def load_env():
    """Load environment variables from project.env file"""
    env_file = Path("project.env")
    
    if env_file.exists():
        with open(env_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#') and '=' in line:
                    key, value = line.split('=', 1)
                    os.environ[key] = value
        print("Environment variables loaded from project.env")
    else:
        print("project.env file not found")

if __name__ == "__main__":
    load_env()
    # Example usage
    git_token = os.environ.get('GIT_AUTH_TOKEN')
    if git_token:
        print(f"Git token loaded: {git_token[:10]}...")
    else:
        print("Git token not found")
