# SSDOT GUI

Designed and built an interactive MATLAB GUI for fNIRS/SSDOT reconstruction, workflow control, and neuroimaging visualization.

## Overview

This project is a MATLAB-based application for guiding users through an SSDOT/fNIRS neuroimaging workflow. It combines a reconstruction pipeline with a custom graphical interface that helps users load data, move through processing steps in order, configure parameters, and visualize reconstruction results.

The goal of this work was to make a complex research workflow more usable, structured, and interactive instead of relying only on scripts and manual processing.

## Problem & Motivation

Neuroimaging pipelines can be difficult to use when they depend on many manual steps, research code, and large data structures. This project improves usability by turning those steps into a guided interface with clearer state management, workflow controls, and visualization tools.

Instead of requiring users to jump between scripts and data files, the GUI helps organize the process from data loading to reconstruction output.

## Features

- Interactive MATLAB GUI for SSDOT/fNIRS workflow navigation
- Step-by-step reconstruction pipeline to help users follow the correct processing order
- Integration with Homer3-derived data and fNIRS file loading
- Data tree organization for browsing runs, sessions, subjects, and group-level outputs
- Shared application state model for storing loaded data, matrices, configuration, and results
- User controls for reconstruction settings such as thresholds, regularization, and device selection
- Visualization tools for sensitivity, spatial basis locations, optical density data, and HbO/HbR results
- Logging and status windows to improve transparency during processing
- Support for averaged result viewing across selected runs or higher-level groupings

## My Contribution

My contribution was the software engineering and GUI development for this project.

I designed and implemented the MATLAB application layer built around the reconstruction pipeline. This included building the interface, structuring the workflow, creating controller logic, developing a class-based shared data model, connecting processing functions to the GUI, and designing tools to prevent users from running steps out of order.

Key areas I worked on include:

- GUI architecture and App Designer integration
- Workflow design for ordered pipeline execution
- Controller logic and event handling
- Class-based shared state management
- Integration of reconstruction functions into the GUI
- Homer3-related data loading and data-tree presentation
- Visualization and user feedback components
- Usability improvements for research-oriented code

## Collaboration

This project was developed in a collaborative, cross-disciplinary research setting. I worked within a biomedical research workflow and focused on translating technical and user needs into a structured MATLAB application. My main contribution centered on software design and application development, including interface design, workflow control, state management, visualization tools, and integration of the reconstruction pipeline into the GUI.

## Tech Stack

- MATLAB
- App Designer
- Homer3 / Homer-derived processing components
- SNIRF-based fNIRS data workflow
- Custom MATLAB classes and visualization utilities

## Project Structure

- `GUI/`  
  Contains the custom GUI, controller logic, shared state model, pipeline windows, logging tools, and visualization helpers.

- `code/`  
  Contains the core reconstruction and processing functions used by the GUI.

- `homer_function/`  
  Contains bundled Homer/SNIRF-related utilities and support code.

- `ST001/`  
  Contains example data, forward-model files, and reconstruction-related assets.

## Skills Demonstrated

- Building research software for end users
- Turning complex technical workflows into guided interfaces
- Connecting frontend-style GUI design with backend processing logic
- Managing shared application state in a multi-step workflow
- Visualizing scientific data in an interactive format
- Contributing effectively in a collaborative technical environment

## Current Status

This project is an active research software effort, and some components are still being refined as the application evolves.

## License

MIT License
