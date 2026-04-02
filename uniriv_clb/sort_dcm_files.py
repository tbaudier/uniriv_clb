#!/usr/bin/env python

import glob
import os
import shutil
from datetime import datetime

import click
import pydicom
import tqdm


def get_tag(ds, tag):
    if tag in ds:
        return ds[tag].value
    else:
        return None


def sort_acquisition(data_dict, nuclear_modality):
    data_dict[nuclear_modality]["acquisition_date"] = sorted(
        data_dict[nuclear_modality]["acquisition_date"], key=lambda x: x[0]
    )
    acquisition_number = 1
    for date, study_instance_uid in data_dict[nuclear_modality]["acquisition_date"]:
        data_dict[nuclear_modality][study_instance_uid] = {
            "cycle": nuclear_modality,
            "tp": "tp" + str(acquisition_number),
        }
        acquisition_number += 1


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--input_folder",
    "-i",
    default=".",
    help="Folder containing DICOM files to sort",
    type=click.Path(exists=True, file_okay=False),
)
@click.option(
    "--output_folder",
    "-o",
    default="./sorted_dcm",
    help="Folder where sorted DICOM files will be stored",
    type=click.Path(file_okay=False),
)
def sort_dicom_files(input_folder, output_folder):
    """
    \b
    Sort DICOM files into a structured directory based on their metadata.
    The directory structure will be:
    output_folder/
        PatientID/
            Modality/
                StudyDescription/
                    StudyDate/
                        StudyTime/
                            SliceThickness/
                                SeriesDescription/
    How to call it:
    python sort_dicom_files.py -i <path to folder containing dcm files> -o <path to output folder where sorted files will be stored>
    """

    # Find all DICOM files in the input folder and its subfolders
    dcm_files = glob.glob(input_folder + "/**/*.[dD][cC][mM]", recursive=True)

    # Initialize the dict containing the count of files per patient and modality
    count_dict = {}

    # Get all acquisition date per modality to classify the files according examen type
    for dcm_file in tqdm.tqdm(dcm_files, desc="Read DICOM files"):
        ds = pydicom.dcmread(dcm_file, stop_before_pixels=True, force=True)
        patient_id = get_tag(ds, (0x0010, 0x0020))  # Patient ID
        if patient_id not in count_dict:
            count_dict[patient_id] = {
                "lu_psma": {
                    "acquisition_date": [],
                },
                "ga": {
                    "acquisition_date": [],
                },
                "fdg": {
                    "acquisition_date": [],
                },
                "other": {
                    "acquisition_date": [],
                },
            }
        modality = get_tag(ds, (0x0008, 0x0060))  # Modality
        if modality == "NM":
            study_instance_uid = get_tag(ds, (0x0020, 0x000D))  # Study Instance UID
            if (
                study_instance_uid not in count_dict[patient_id]
                or count_dict[patient_id][study_instance_uid] == "other"
            ):
                count_dict[patient_id][study_instance_uid] = "lu_psma"
                acquisition_date = get_tag(ds, (0x0008, 0x0022))  # Acquisition Date
                count_dict[patient_id]["lu_psma"]["acquisition_date"].append(
                    (acquisition_date, study_instance_uid)
                )
        elif modality == "PT":
            radiopharmaceutical = ds[0x0054, 0x0016][0][0x0054, 0x0300][0][
                0x0008, 0x0104
            ].value
            study_instance_uid = get_tag(ds, (0x0020, 0x000D))  # Study Instance UID
            if "68" in radiopharmaceutical.lower():
                if (
                    study_instance_uid not in count_dict[patient_id]
                    or count_dict[patient_id][study_instance_uid] == "other"
                ):
                    count_dict[patient_id][study_instance_uid] = "ga"
                    acquisition_date = get_tag(ds, (0x0008, 0x0022))  # Acquisition Date
                    count_dict[patient_id]["ga"]["acquisition_date"].append(
                        (acquisition_date, study_instance_uid)
                    )
            elif "18" in radiopharmaceutical.lower():
                if (
                    study_instance_uid not in count_dict[patient_id]
                    or count_dict[patient_id][study_instance_uid] == "other"
                ):
                    count_dict[patient_id][study_instance_uid] = "fdg"
                    acquisition_date = get_tag(ds, (0x0008, 0x0022))  # Acquisition Date
                    count_dict[patient_id]["fdg"]["acquisition_date"].append(
                        (acquisition_date, study_instance_uid)
                    )
            else:
                if study_instance_uid not in count_dict[patient_id]:
                    count_dict[patient_id][study_instance_uid] = "other"
        else:
            study_instance_uid = get_tag(ds, (0x0020, 0x000D))  # Study Instance UID
            if study_instance_uid not in count_dict[patient_id]:
                count_dict[patient_id][study_instance_uid] = "other"
    for patient_id in count_dict:
        for siu in count_dict[patient_id]:
            if count_dict[patient_id][siu] == "other":
                acquisition_date = get_tag(ds, (0x0008, 0x0022))  # Acquisition Date
                count_dict[patient_id]["other"]["acquisition_date"].append(
                    ("20200101", siu)
                )

    # Sort the acquisition date per modality
    for patient_id in count_dict:
        sort_acquisition(count_dict[patient_id], "fdg")
        sort_acquisition(count_dict[patient_id], "ga")
        sort_acquisition(count_dict[patient_id], "other")
        count_dict[patient_id]["lu_psma"]["acquisition_date"] = sorted(
            count_dict[patient_id]["lu_psma"]["acquisition_date"], key=lambda x: x[0]
        )
        cycle_number = 1
        acquisition_number = 1
        previous_date = None
        for date, study_instance_uid in count_dict[patient_id]["lu_psma"][
            "acquisition_date"
        ]:
            if previous_date is None:
                acquisition_number = 1
                count_dict[patient_id]["lu_psma"][study_instance_uid] = {
                    "cycle": "cycle" + str(cycle_number),
                    "tp": "tp" + str(acquisition_number),
                }
                cycle_number += 1
                acquisition_number += 1
                previous_date = date
            elif (
                datetime.strptime(date, "%Y%m%d").date()
                - datetime.strptime(previous_date, "%Y%m%d").date()
            ).days > 30:
                acquisition_number = 1
                count_dict[patient_id]["lu_psma"][study_instance_uid] = {
                    "cycle": "cycle" + str(cycle_number),
                    "tp": "tp" + str(acquisition_number),
                }
                cycle_number += 1
                acquisition_number += 1
                previous_date = date
            else:
                count_dict[patient_id]["lu_psma"][study_instance_uid] = {
                    "cycle": "cycle" + str(cycle_number - 1),
                    "tp": "tp" + str(acquisition_number),
                }
                acquisition_number += 1

    # Copy dcm files into the appropriate folder based on timepoint and cycle for each modality
    for dcm_file in tqdm.tqdm(dcm_files, desc="Classifying DICOM files"):
        ds = pydicom.dcmread(dcm_file, stop_before_pixels=True, force=True)
        study_instance_uid = get_tag(ds, (0x0020, 0x000D))  # Study Instance UID
        patient_id = get_tag(ds, (0x0010, 0x0020))  # Patient ID
        if study_instance_uid in count_dict[patient_id]:
            folder_dict = count_dict[patient_id][
                count_dict[patient_id][study_instance_uid]
            ][study_instance_uid]
            modality = get_tag(ds, (0x0008, 0x0060))  # Modality
            if modality == "CT":
                modality_folder = "ct"
            elif modality == "NM" or modality == "PT":
                modality_folder = "nm"
            new_dir = os.path.join(
                str(output_folder),
                str(patient_id),
                str(folder_dict["cycle"]),
                str(folder_dict["tp"]),
                str(modality_folder),
                "dcm",
            )
            os.makedirs(new_dir, exist_ok=True)
            filename = os.path.basename(dcm_file)
            new_path = os.path.join(new_dir, filename)
            shutil.copy(dcm_file, new_path)


if __name__ == "__main__":
    sort_dicom_files()
