# Revision Workflow Reminder

## Critical: Dual Repository Workflow

When making changes to address reviewer comments, **ALWAYS** work with TWO repositories:

### 1. overleaf-paper repository (main paper files)
- **Location**: `/data/yizhangh/ldp-pq/overleaf-paper/`
- **Purpose**: Contains the actual paper files (1intro.tex, 2prelim.tex, etc.)
- **Remote**: Overleaf paper repository

### 2. overleaf-revision-letter repository (revision responses)  
- **Location**: `/data/yizhangh/ldp-pq/overleaf-revision-letter/`
- **Purpose**: Contains reviewer response files (reviewer4.tex, etc.)
- **Remote**: Overleaf revision letter repository

## Required Steps for Each Change

For every revision, follow this complete workflow:

### Step 1: Make Changes to Main Paper
```bash
cd /data/yizhangh/ldp-pq/overleaf-paper/
# Edit files (e.g., 1intro.tex)
git add <modified-files>
git commit -m "Descriptive commit message"
git push
```

### Step 2: Update Revision Letter Response
```bash
cd /data/yizhangh/ldp-pq/overleaf-revision-letter/
# Update corresponding response files (e.g., reviewer4.tex)
git add <modified-files>
git commit -m "Update reviewer response to reflect implemented changes"
git push
```

## Example Workflow

**Problem**: [Description of reviewer comment]

**Solution Applied**:
1. **overleaf-paper/[relevant-file].tex**: [Description of changes made to main paper]
2. **overleaf-revision-letter/[relevant-file].tex**: [Description of response updates]

## Key Reminders

- **Both repositories must be updated** for each revision
- **Both repositories must be committed and pushed** to their respective remotes
- **Revision letter responses should reflect completed changes**, not just promises
- **Use descriptive commit messages** explaining what was changed and why
- **Highlight modified content in blue** for easy tracking
- **In revision letters, make all section/figure/table references blue** using `{\color{blue}Section X}`, `{\color{blue}Figure~Y}`, `{\color{blue}Table Z}`, etc.

## AI Assistant Guidelines

### When User Says "Address Comment X"
1. **Read the specific reviewer comment** from the revision letter file
2. **Identify which files need updating**:
   - Main paper files in `/data/yizhangh/ldp-pq/overleaf-paper/`
   - Revision response in `/data/yizhangh/ldp-pq/overleaf-revision-letter/`
3. **Make changes to main paper first** with blue coloring
4. **Update revision letter response** to reflect completed work (not promises)
5. **Commit and push both repositories** automatically

### File Location Quick Reference
- **Main paper files**: `/data/yizhangh/ldp-pq/overleaf-paper/[filename].tex`
- **Revision responses**: `/data/yizhangh/ldp-pq/overleaf-revision-letter/reviewer[X].tex`
- **Main paper repo**: `/data/yizhangh/ldp-pq/overleaf-paper/` (independent git repo)
- **Revision letter repo**: `/data/yizhangh/ldp-pq/overleaf-revision-letter/` (independent git repo)

### Always Do These Steps
1. **Read the reviewer comment** to understand the issue
2. **Find relevant files** in main paper that need updates
3. **Make changes with blue coloring** for new/modified content
4. **Update revision letter response** to reflect completed work
5. **Commit and push both repositories** with descriptive messages
6. **Never leave responses as "will do"** - always show completed work

### Critical Rule for Revision Letters
- **NEVER say "we will xxx"** - always say "we have done xxx"
- **ALWAYS specify where**: "we have clarified in Section X" or "we have updated Figure Y"
- **Use blue coloring for all references**: `{\color{blue}Section X}`, `{\color{blue}Figure~Y}`, `{\color{blue}Table Z}`
- **DO NOT color entire response text** - only color section/figure/table references
- **Be specific about what was actually revised**: When writing responses, consider what we actually revised in the main paper and specify exactly how and where the changes were made

**Note**: This document has evolved from a revision plan (where "we will do X" was appropriate) to a revision letter (where "we have done X" is required). Revision letters must reflect completed work, not future promises.

## Response Writing Guidelines

### Tone and Style
- **Keep responses concise and direct** - avoid overly apologetic language
- **Match the user's professional tone** - not verbose or wordy
- **Use "we have clarified/changed" instead of "we will clarify/change"**
- **Place key actions at the beginning** of responses for clarity

### Structure for Responses
1. **Thank the reviewer** (brief)
2. **State what was done** (with blue section references)
3. **Provide brief explanation** of the changes
4. **Reference prior work** when relevant (for justification)
5. **Be specific about revisions**: Describe exactly what changes were made in the main paper, not just "we clarified" - specify the mechanism, purpose, and location

### Blue Coloring Usage
- **All new/modified content in main paper**: `{\color{blue}...}`
- **Section references in revision letters**: `{\color{blue}Section X}`
- **Figure/table references**: `{\color{blue}Figure~Y}`, `{\color{blue}Table Z}`
- **For revision letter responses**: Only color section/figure/table references, NOT the entire response text

## Repository Management Notes

### Important: Main Repository vs Sub-repositories
- **Main repository** (`/data/yizhangh/ldp-pq/`) does NOT track the overleaf sub-repositories
- **Each overleaf repository is independent** and must be managed separately
- **No need to update main repository** when working on paper/revision changes
- **Always work directly in the specific overleaf repository** for that type of content
- **NEVER commit overleaf repos from main directory**: Never run `cd /data/yizhangh/ldp-pq/overleaf-* && git push` - always work directly in the overleaf repository directories

### Repository URLs
- **Main paper**: `https://git.overleaf.com/67bffc24e87bec767040ba16`
- **Revision letter**: `https://git.overleaf.com/68b7b08f2c7f5dd8ff0da856`

## Common Mistakes to Avoid

❌ **Only updating one repository**
❌ **Forgetting to push changes to remote**
❌ **Leaving revision responses as promises instead of completed actions**
❌ **Not documenting what specific changes were made**
❌ **Not making section/figure/table references blue in revision letters**
❌ **Trying to track overleaf repos from main repository**

✅ **Always update both repositories**
✅ **Always commit and push both repositories**
✅ **Always update revision responses to reflect completed work**
✅ **Always use descriptive commit messages**
✅ **Always make section/figure/table references blue in revision letters**
✅ **Always work directly in the specific overleaf repository**
