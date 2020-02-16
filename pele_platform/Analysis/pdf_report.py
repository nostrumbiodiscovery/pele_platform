from fpdf import FPDF
import glob
import os


def create_report(plots, clusters, top_poses, output="simulation_summary.pdf"):


    pdf = FPDF()
    pdf.add_page() 

    # Title
    pdf.set_font('Arial', 'B', 15)
    pdf.cell(0, 10, "Score {} kcal/mol", align='C')

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
    for i, poses in enumerate(top_poses[0:10]):
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

if __name__ == "__main__":
    plots = glob.glob("/work/NBD_Utilities/PELE/PELE_Softwares/pele_platform_devel/pele_platform/results/Plots/*")
    clusters = glob.glob("/work/NBD_Utilities/PELE/PELE_Softwares/pele_platform_devel/pele_platform/results/clusters/*.png")
    top_poses = glob.glob("/work/NBD_Utilities/PELE/PELE_Softwares/pele_platform_devel/pele_platform/results/BestStructs/*.pdb")
    create_report(plots, clusters, top_poses)
