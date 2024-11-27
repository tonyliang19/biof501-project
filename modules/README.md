# Modules files

This directory stores nextflow module codes that stores process definitions like input, output, and other process metadata. This way, we could decouple process logic from the main workflow logics, leaving all relevant implementations here.

There two subdirectories here: `local` and `nf-core`. The first one are locally created custom modules used in the pipeline, and the second one being a generic utility module widely used for nf-core pipelines, specifically for saving software versions, and logging purposes.

The naming conventions and directory structure of `local` modules are:

```bash
# tool_name is the binary used for computation
# main.nf is the process definition of the tool, specifying what input/output it contains,
# actual script that is being run and other metadata like where to publish result, container used, etc. 
<tool_name>
  └── main.nf
# Some tools that could be break down into subtools might have this naming:
<big_tool_name>/
├── <big_tool-sub_tool1>
│   └── main.nf
└── <big_tool-sub_tool2>
    └── main.nf
```