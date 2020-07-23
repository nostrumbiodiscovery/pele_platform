from fpdf import FPDF
import os
import numpy as np


def create_report(plots, clusters, top_poses, best_energies, output="simulation_summary.pdf"):


    pdf = FPDF()
    pdf.add_page() 

    # Title
    pdf.set_font('Arial', 'B', 15)
    idx_25 = int(len(best_energies)*0.25)+1
    mean_25_quartile = np.mean(sorted(best_energies)[:idx_25])
    pdf.cell(0, 10, "Score {} kcal/mol".format(mean_25_quartile), align='C')

    pdf.ln(40) #4 linebreaks


    #Plots
    pdf.set_font('Arial', 'B', size=12)
    pdf.cell(0, 10, "Plots", align='C')
    pdf.set_font('Arial', size=10)
    for i, plot in enumerate(plots):
        if i%2 == 0 and i!=0:
            pdf.ln(1000000) #pagebreak
        pdf.ln(10)
        name = os.path.basename(plot)[:-4]
        pdf.cell(0, 50, name, align='C')
        pdf.ln(30)
        pdf.image(plot, x=70, w=83)

    pdf.ln(10000000) #page break


    #Top poses
    pdf.set_font('Arial', 'B', size=12)
    pdf.cell(0, 10, "Top poses", align='C')
    pdf.set_font('Arial', size=10)
    top_poses_ordered = np.array(top_poses)[np.argsort(best_energies)[:len(top_poses)]]
    for i, poses in enumerate(top_poses_ordered[0:20]):
        if i == 0:
            pdf.ln(10)
        pdf.cell(0, 50, os.path.basename(poses), align='C')
        pdf.ln(10)

    pdf.ln(10000000) #page break

    #Clusters
    pdf.set_font('Arial', 'B', size=12)
    pdf.cell(0, 10, "Clusters", align='C')
    pdf.set_font('Arial', size=10)
    for i, plot in enumerate(clusters):
        if i%2 == 0 and i!=0:
            pdf.ln(60)
        pdf.ln(10)
        name = os.path.basename(plot)[:-4]
        pdf.cell(0, 50, name, align='C')
        pdf.ln(30)
        pdf.image(plot, x=70, w=83)

    #Output report    
    pdf.output(output, 'F')
    return output
