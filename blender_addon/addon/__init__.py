bl_info = {
    "name": "SPH Solver",
    "author": "José Manuel Valverde Pérez",
    "version": (1, 0),
    "blender": (2, 93, 0),
    "location": "Properties > Physics",
    "description": "SPH solver for Blender scenes",
    "warning": "",
    "wiki_url": "",
    "tracker_url": "",
    "category": "Physics"
}

import bpy
from bpy.props import StringProperty, IntProperty, CollectionProperty
from bpy.types import PropertyGroup, UIList, Operator, Panel

mode_options = [
    ("wcsph_method", "Weakly Compressible SPH", '', '', 0),
    ("pcisph_method", "Predictive-Corrective Incompressible SPH", '', '', 1),
    ("dfsph_method", "Divergence-Free SPH", '', '', 2)
]

fluid_unit_types = [
    ("fu_fluid_block_type", "Fluid block", '', '', 0),
    ("fu_geometry_type", "Geometry", '', '', 2),
    ("fu_emitter_type", "Emitter", '', '', 1)
]

boundary_unit_types = [
    ("bu_box_type", "Box", '', '', 0),
    ("bu_sphere_type", "Sphere", '', '', 2),
    ("bu_geometry_type", "Geometry", '', '', 1)
]

boundary_options = [
    ("none_boundary_method", "None", '', '', 0),
    ("cube_boundary_method", "Simple boundary handling", '', '', 1),
    ("pcisph_boundary_method", "PCISPH boundary handling", '', '', 2),
    ("akinci_boundary_method", "Akinci boundary handling", '', '', 3)
]

emitter_types = [
    ("none_emitter", "None", '', '', 0),
    ("square_emitter", "Square emitter", '', '', 1),
    ("circle_emitter", "Circle emitter", '', '', 2)
]

viscosity_options = [
    ("none_visco_method", "None", '', '', 0),
    ("standard_visco_method", "Standard Viscosity", '', '', 1),
    ("artificial_visco_method", "Artificial Viscosity", '', '', 2),
    ("xsph_visco_method", "XSPH Viscosity", '', '', 3)
]

surften_options = [
    ("none_st_method", "None", '', '', 0),
    ("akinci_st_method", "Akinci Surface Tension", '', '', 1)
]

adhesion_options = [
    ("none_adhesion_method", "None", '', '', 0),
    ("akinci_adhesion_method", "Akinci Adhesion", '', '', 1)
]

######## JSON IO #################

import json
import os
from bpy_extras.io_utils import ImportHelper

def read_frames(context, files, path):
    
    path = os.path.split(path)[0]
    filepath = os.path.join(path, files[0].name)
    numFiles = len(files)
    
    # Leer el primer frame para saber el numero de particulas
    numParticles = 0
    numActiveParticles = 0
    try:
        with open(filepath, 'r') as infile:
            state = json.load(infile)
            numParticles = state["numParticles"] 
            numActiveParticles = state["numActiveParticles"]
    except Exception as error:
        print("Cannot start read", filepath)
    
    if numParticles > 0:
        context.scene.frame_end = len(files)

        ps_index = -1
        keys = context.object.particle_systems.keys()
        for i in range(0, len(keys)):
            if "SPHSystem" == keys[i]:
                ps_index = i
                break
            
        if ps_index > -1:
            p_sist = context.object.particle_systems['SPHSystem']
            context.object.particle_systems.active_index = ps_index
            bpy.ops.object.particle_system_remove()
            
        p_sist = context.object.modifiers.new("SPHSystem", 'PARTICLE_SYSTEM').particle_system
            
        # Establecer settings del sistema de particulas
        p_sist.settings.count = numParticles
        p_sist.settings.physics_type = 'NEWTON'
        p_sist.settings.emit_from = 'VERT'
        p_sist.settings.timestep = 0
        p_sist.settings.mass = 0
        p_sist.settings.normal_factor = 0
        p_sist.settings.frame_end = p_sist.settings.frame_start = 1
        p_sist.settings.lifetime = len(files)
        p_sist.particles.data.settings.display_size = context.object.sph_scene_data.particleRadius
        p_sist.settings.display_color = 'VELOCITY'
        p_sist.settings.color_maximum = 4
        
        # Limpiar cache
        fake_context = context.copy()
        fake_context["point_cache"] = p_sist.point_cache
        bpy.ops.ptcache.free_bake(fake_context)
 
        # Necesario para coger el psist actualizado
        degp = context.evaluated_depsgraph_get()
        p_sist = context.object.evaluated_get(degp).particle_systems[p_sist.name]
        
        # velocidades a 0 para que el integrador de blender no haga nada
        context.scene.frame_set(0)
        vel = [0,0,0] * numParticles
        p_sist.particles.foreach_set('velocity',vel)

        print("Start reading frames...")
        for file in files:
            filepath = os.path.join(path, file.name)
            split = file.name.split("_")
            split = split[1].split(".")
            frame = int(split[0])
            try:
                with open(filepath, 'r') as infile:
                    state = json.load(infile)
                    print("Inserting frame", str(frame) + "/" + str(numFiles), "in particle system", end="\r")
                    insert_frame(context.scene, frame, state, p_sist) 
            except Exception as error:
                stop = True
                print("")
                print("Cannot read", filepath)
        
        print("")  
        print("Frames read!")
        
        # bake cache
        bpy.ops.ptcache.bake_from_cache(fake_context)
        context.scene.frame_set(1)                
            

def insert_frame(scene, frame, state, p_sist): 
    # Numap debe ser menor que nump
    nump = state["numParticles"]
    numap = state["numActiveParticles"]
    
    scene.frame_set(frame)
    
    locations = [0,0,0] * nump
    for i in range(0, 3 * numap, 3):
        locations[i] = state["positions"][i//3][0]
        locations[i+1] = state["positions"][i//3][1]
        locations[i+2] = state["positions"][i//3][2]
        
    velocities = [0,0,0] * nump
    for i in range(0, 3 * numap, 3):
        velocities[i] = state["velocities"][i//3][0]
        velocities[i+1] = state["velocities"][i//3][1]
        velocities[i+2] = state["velocities"][i//3][2]
        
    # Hide inactive particles
    for i in range(3 * numap, 3 * nump, 3):
        locations[i] = 1000
        locations[i+1] = 1000
        locations[i+2] = 1000
    
    p_sist.particles.foreach_set("location", locations)
    p_sist.particles.foreach_set("velocity", velocities)
    

def import_config(object, filepath):

    with open(filepath, 'r') as infile:
        config = json.load(infile)
    
        print("Import scene data")
        import_scene(object, config)
        
        print("Import fluid set data")
        import_fluid_set(object, config)
        
        print("Import boundary set")
        import_boundary_set(object, config)
        
    return

def import_scene(object, config):
    
    sceneData = object.sph_scene_data
    boundarySet = object.sph_boundary_set
    
    simulationMethod = config["Configuration"]["simulationMethod"]        
    sceneData.simulationMethod = mode_options[simulationMethod][0]

    boundaryMethod = config["Configuration"]["boundaryMethod"]
    boundarySet.boundaryMethod = boundary_options[boundaryMethod + 1][0]
    
    sceneData.startTime = config["Configuration"]["startTime"]
    sceneData.endTime = config["Configuration"]["endTime"]
    sceneData.timeStep = config["Configuration"]["timeStep"]
    sceneData.fps = config["Configuration"]["fps"]
    sceneData.minTimeStep = config["Configuration"]["minTimeStep"]
    sceneData.maxTimeStep = config["Configuration"]["maxTimeStep"]
    sceneData.particleRadius = config["Configuration"]["particleRadius"]
    sceneData.gravity = config["Configuration"]["gravity"]
    
    if sceneData.simulationMethod == mode_options[0][0]: # WCSPH
        sceneData.stiffness = config['Configuration']['stiffness']
        sceneData.gamma = config['Configuration']['gamma']
        sceneData.cflFactor = config['Configuration']['cflFactor']
    elif sceneData.simulationMethod == mode_options[1][0]: # PCISPH
        sceneData.etaPCI = config['Configuration']['eta']
        sceneData.minIterationsPCI = config['Configuration']['minIterations']
        sceneData.maxIterationsPCI = config['Configuration']['maxIterations']
        sceneData.cflFactor = config['Configuration']['cflFactor']
    elif sceneData.simulationMethod == mode_options[2][0]: # DFSPH
        sceneData.eta = config['Configuration']['eta']
        sceneData.minIterationsDF = config['Configuration']['minIterations']
        sceneData.maxIterationsDF = config['Configuration']['maxIterations']
        sceneData.etaV = config['Configuration']['etaV']
        sceneData.minIterationsDFV = config['Configuration']['minIterations']
        sceneData.maxIterationsDFV = config['Configuration']['maxIterations']
        sceneData.cflFactor = config['Configuration']['cflFactor']
    
def import_fluid_set(object, config):
    
    fluidSet = object.sph_fluid_set
    
    if "Fluid" in config.keys():
        
        viscosityMethod = config["Fluid"]["viscosityMethod"]
        fluidSet.viscosityMethod = viscosity_options[viscosityMethod + 1][0]

        surftenMethod = config["Fluid"]["surfaceTensionMethod"]
        if surftenMethod == -1:
            fluidSet.surftenMethod = surften_options[0][0]
        elif surftenMethod == 3:
            fluidSet.surftenMethod = surften_options[1][0]
        
        adhesionMethod = config["Fluid"]["adhesionMethod"]
        if adhesionMethod == -1:
            fluidSet.adhesionMethod == adhesion_options[0][0]
        elif adhesionMethod == 4:
            fluidSet.adhesionMethod == adhesion_options[1][0]
            
        # Fluid models
        for i in range(0, len(config["Fluid"]["FluidModel"])):
            fm = fluidSet.fs_list.add()
            fm.density0 = config["Fluid"]["FluidModel"][i]["density0"]
            fm.viscosity = config["Fluid"]["FluidModel"][i]["viscosity"]
            fm.bviscosity = config["Fluid"]["FluidModel"][i]["boundaryViscosity"]
            fm.surften = config["Fluid"]["FluidModel"][i]["surfaceTension"]
            fm.adhesion = config["Fluid"]["FluidModel"][i]["adhesion"]
            
            if "fluidBlocks" in config["Fluid"]["FluidModel"][i].keys():
                fluidBlocks = config["Fluid"]["FluidModel"][i]["fluidBlocks"]
                for j in range(0, len(fluidBlocks)):
                    fu = fm.fm_list.add()
                    
                    fu.min = fluidBlocks[j]["min"]
                    fu.max = fluidBlocks[j]["max"]
                    
                    x = (fu.max[0] + fu.min[0]) * 0.5
                    y = (fu.max[1] + fu.min[1]) * 0.5
                    z = (fu.max[2] + fu.min[2]) * 0.5
                    
                    dx = fu.max[0] - fu.min[0]
                    dy = fu.max[1] - fu.min[1]
                    dz = fu.max[2] - fu.min[2]
                    
                    # crear el fluid block
                    # guardar objeto activo antes de crear el nuevo cubo
                    nameActiveObject = bpy.context.view_layer.objects.active.name

                    # crear el cubo
                    bpy.ops.mesh.primitive_cube_add(size=0.5, enter_editmode=False, align="WORLD", location=(0, 0, 0), scale=(1, 1, 1))
                    bpy.context.view_layer.objects.active.name = "Fluid block"
                    fluidBlock = bpy.context.view_layer.objects.active
                    print("FluidBlock importado")
                    
                    # si no se establece el vector de golpe solo se asigna correctamente la ultima componente asignada
                    fluidBlock.location = (x, y, z)
                    fluidBlock.dimensions = (dx, dy, dz)

                    # Guardar en el fluid unit una referencia al objeto que lo representa y su tipo
                    fu.object = bpy.context.view_layer.objects.active
                    fu.fu_type = fluid_unit_types[0][0]

                    # reestablecer el objeto activo
                    bpy.context.view_layer.objects.active = bpy.data.objects[nameActiveObject]
            
            #elif "geometry" in config["Fluid"]["FluidModel"][i].keys():
            
            elif "emitters" in config["Fluid"]["FluidModel"][i].keys():
                emitters = config["Fluid"]["FluidModel"][i]["emitters"]
                
                for j in range(0, len(emitters)):
                    fu = fm.fm_list.add()
                    
                    fu.velocity = emitters[j]["velocity"]
                    fu.numParticles = emitters[j]["numParticles"]
                    fu.startTime = emitters[j]["startTime"]
                    fu.spacing = emitters[j]["spacing"]
                    
                    type = emitters[j]["type"]
                    fu.emitter_type = emitter_types[type + 1][0] # al establecer el tipo es cuando se crean los emitters por lo que no hay que crearlos manualmente
                    
                    if type == 0: # square
                        fu.object.dimensions = (emitters[j]["width"], emitters[j]["height"], 0)
                        fu.empty.scale.z = 0.5 * (emitters[j]["width"] + emitters[j]["height"])
                        
                    elif type == 1: # circle
                        fu.object.dimensions = (emitters[j]["width"] * 2.0, emitters[j]["width"] * 2.0, 0)
                        fu.empty.scale.z = emitters[j]["width"] * 2.0
                        
                    fu.object.location = (emitters[j]["position"][0], emitters[j]["position"][1], emitters[j]["position"][2])
                    
                    rotMode = fu.object.rotation_mode
                    fu.object.rotation_mode = 'QUATERNION'
                    fu.object.rotation_quaternion = (emitters[j]["rotation"][0], emitters[j]["rotation"][1], emitters[j]["rotation"][2], emitters[j]["rotation"][3])
                    fu.object.rotation_mode = rotMode
    else:
        # si no hay fluid, crear un fluid model ya que si no el addon se iniciara en no multifase sin ningun fm
        fluidSet = object.sph_fluid_set.fs_list.add()
        
def import_boundary_set(object, config):
    
    boundarySet = object.sph_boundary_set
    
    if "Boundary" in config.keys():
        
        # Boundary models
        for i in range(0, len(config["Boundary"])):
            bm = boundarySet.bs_list.add()
            
            if "box" in config["Boundary"][i].keys():
                boxes = config["Boundary"][i]["box"]
                for j in range(0, len(boxes)):
                    bu = bm.bm_list.add()
                    
                    bu.min = boxes[j]["min"]
                    bu.max = boxes[j]["max"]
                    
                    x = (bu.max[0] + bu.min[0]) * 0.5
                    y = (bu.max[1] + bu.min[1]) * 0.5
                    z = (bu.max[2] + bu.min[2]) * 0.5
                    
                    dx = bu.max[0] - bu.min[0]
                    dy = bu.max[1] - bu.min[1]
                    dz = bu.max[2] - bu.min[2]
      
                    # crear el box
                    # guardar objeto activo antes de crear el nuevo cubo
                    nameActiveObject = bpy.context.view_layer.objects.active.name

                    # crear el cubo
                    bpy.ops.mesh.primitive_cube_add(size=0.5, enter_editmode=False, align="WORLD", location=(0, 0, 0), scale=(1, 1, 1))
                    bpy.context.view_layer.objects.active.name = "Boundary cube"
                    box = bpy.context.view_layer.objects.active
                    print("box importado")
                    
                    # si no se establece el vector de golpe solo se asigna correctamente la ultima componente asignada
                    box.location = (x, y, z)
                    box.dimensions = (dx, dy, dz)

                    # Guardar en el fluid unit una referencia al objeto que lo representa y su tipo
                    bu.object = bpy.context.view_layer.objects.active
                    bu.bu_type = boundary_unit_types[0][0]

                    # reestablecer el objeto activo
                    bpy.context.view_layer.objects.active = bpy.data.objects[nameActiveObject]
                    
            if "sphere" in config["Boundary"][i].keys():
                spheres = config["Boundary"][i]["sphere"]
                for j in range(0, len(spheres)):
                    bu = bm.bm_list.add()
                    
                    # crear el sphere
                    # guardar objeto activo antes de crear el nuevo cubo
                    nameActiveObject = bpy.context.view_layer.objects.active.name

                    # crear el cubo
                    bpy.ops.mesh.primitive_uv_sphere_add(radius=0.5, enter_editmode=False, align="WORLD", location=(0, 0, 0), scale=(1, 1, 1))
                    bpy.context.view_layer.objects.active.name = "Boundary sphere"
                    sphere = bpy.context.view_layer.objects.active
                    print("sphere importado")
                    
                    # si no se establece el vector de golpe solo se asigna correctamente la ultima componente asignada
                    location = spheres[j]["pos"]
                    radius = spheres[j]["radius"]
                    sphere.location = (location[0], location[1], location[2])
                    sphere.dimensions = (radius * 2.0, radius * 2.0, radius * 2.0)

                    # Guardar en el fluid unit una referencia al objeto que lo representa y su tipo
                    bu.object = bpy.context.view_layer.objects.active
                    bu.bu_type = boundary_unit_types[1][0]

                    # reestablecer el objeto activo
                    bpy.context.view_layer.objects.active = bpy.data.objects[nameActiveObject]
                    
            if "geometry" in config["Boundary"][i].keys():
                geometries = config["Boundary"][i]["geometry"]
                
                for j in range(0, len(geometries)):
                    
                    spacing = geometries[j]["spacing"]
                    path = geometries[j]["path"]
                    
                    try:
                        bpy.ops.import_scene.obj(filepath=path)
                    except (RuntimeError, FileNotFoundError) as e:
                        print("ERROR: File " + path + " not found")
                        dir = os.path.split(path)[0]
                        line1 = "File " + path + " not found."
                        line2 =  "Put the .obj in the directory " + dir + " or add it manually."
                        showMessageBox(title = "File not found", icon = 'ERROR', lines = (line1, line2))
                    
                    #coger el objeto importado y asignarlo a bu.object
                    bu = bm.bm_list.add()    
                    bu.bu_type = boundary_unit_types[2][0]
                

def export_config(object):

    config = {}

    # Scene data
    print("Export scene data")
    export_scene(object, config)

    # Fluid set
    print("Export fluid set")
    export_fluid_set(object, config)

    # Boundary set
    print("Export boundary set")
    export_boundary_set(object, config)
    
    filename = os.path.split(bpy.data.filepath)[1]
    filename = filename[:-6] + ".json"

    with open(filename, 'w') as outfile:
        json.dump(config, outfile, indent=4)
        print("Exported scene")
        outfile.close()

    return {"FINISHED"}

def export_scene(object, config):

    sceneData = object.sph_scene_data
    boundarySetData = object.sph_boundary_set

    simulationMethod = -1
    if sceneData.simulationMethod == mode_options[0][0]: # WCSPH
        simulationMethod = 0
    elif sceneData.simulationMethod == mode_options[1][0]: # PCISPH
        simulationMethod = 1
    elif sceneData.simulationMethod == mode_options[2][0]: # DFSPH
        simulationMethod = 2

    boundaryMethod = -1
    if boundarySetData.boundaryMethod == boundary_options[1][0]: # Simple cube bh
        boundaryMethod = 0
    elif boundarySetData.boundaryMethod == boundary_options[2][0]: # PCISPH bh
        boundaryMethod = 1
    if boundarySetData.boundaryMethod == boundary_options[3][0]: # Akinci bh
        boundaryMethod = 2

    config['Configuration'] = {

        "startTime" : sceneData.startTime,
        "endTime":sceneData.endTime,
        "timeStep":sceneData.timeStep,
        "fps":sceneData.fps,
        "minTimeStep":sceneData.minTimeStep,
        "maxTimeStep":sceneData.maxTimeStep,

        "gravity":[sceneData.gravity[0],sceneData.gravity[1],sceneData.gravity[2]],
        "particleRadius":sceneData.particleRadius,
        "simulationMethod":simulationMethod,
        "boundaryMethod":boundaryMethod,
    }

    # Solver data
    if sceneData.simulationMethod == mode_options[0][0]: # WCSPH
        config['Configuration']['stiffness'] = sceneData.stiffness
        config['Configuration']['gamma'] = sceneData.gamma
        config['Configuration']['cflFactor'] = sceneData.cflFactor
    elif sceneData.simulationMethod == mode_options[1][0]: # PCISPH
        config['Configuration']['eta'] = sceneData.etaPCI
        config['Configuration']['minIterations'] = sceneData.minIterationsPCI
        config['Configuration']['maxIterations'] = sceneData.maxIterationsPCI
        config['Configuration']['cflFactor'] = sceneData.cflFactor
    elif sceneData.simulationMethod == mode_options[2][0]: # DFSPH
        config['Configuration']['eta'] = sceneData.eta
        config['Configuration']['minIterations'] = sceneData.minIterationsDF
        config['Configuration']['maxIterations'] = sceneData.maxIterationsDF
        config['Configuration']['etaV'] = sceneData.etaV
        config['Configuration']['minIterations'] = sceneData.minIterationsDFV
        config['Configuration']['maxIterations'] = sceneData.maxIterationsDFV
        config['Configuration']['cflFactor'] = sceneData.cflFactor

def export_fluid_set(object, config):

    fluidSetData = object.sph_fluid_set

    viscosityMethod = -1
    if fluidSetData.viscosityMethod == viscosity_options[1][0]: # Standard visco
        viscosityMethod = 0
    elif fluidSetData.viscosityMethod == viscosity_options[2][0]: # Artificial visco
        viscosityMethod = 1
    elif fluidSetData.viscosityMethod == viscosity_options[3][0]: # XSPH visco
        viscosityMethod = 2

    surftenMethod = -1
    if fluidSetData.surftenMethod == surften_options[1][0]: # Akinci surface tension
        surftenMethod = 3

    adhesionMethod = -1
    if fluidSetData.adhesionMethod == adhesion_options[1][0]: # Adhesion
        adhesionMethod = 4

    fluidSet = {}

    for fm in fluidSetData.fs_list:

        if len(fm.fm_list) > 0: 

            fluidSet['FluidModel'] = []
                
            fluidBlocks = []
            #geometry = []
            emitters = []

            for fu in fm.fm_list:

                if fu.fu_type == fluid_unit_types[0][0]: # fluid block

                    sideX = fu.object.dimensions[0]
                    sideY = fu.object.dimensions[1]
                    sideZ = fu.object.dimensions[2]

                    x = fu.object.location[0]
                    y = fu.object.location[1]
                    z = fu.object.location[2]

                    min = [x - sideX * 0.5,
                           y - sideY * 0.5,
                           z - sideZ * 0.5]
                    max = [x + sideX * 0.5,
                           y + sideY * 0.5,
                           z + sideZ * 0.5]

                    fluidBlock = {
                        "min":min,
                        "max":max
                    }

                    fluidBlocks.append(fluidBlock)

                #elif fu.fu_type == fluid_unit_types[1][0]: # geometry

                    # Hay que exportar las geometrias en formato obj en determinada ruta
                    # Aun no se ha implementado el sampleo de volumen en el simulador por lo que de momento nada

                elif fu.fu_type == fluid_unit_types[2][0]: # emitter

                    if fu.emitter_type != emitter_types[0][0]: # si es distinto de none

                        emitter = {
                            "velocity":fu.velocity,
                            "numParticles":fu.numParticles,
                            "startTime":fu.startTime,
                            "spacing":fu.spacing
                        }

                        x = fu.object.location[0]
                        y = fu.object.location[1]
                        z = fu.object.location[2]

                        emitter['position'] = [x, y, z]
                        
                        quat = fu.object.rotation_euler.to_quaternion()

                        w = quat[0]
                        x = quat[1]
                        y = quat[2]
                        z = quat[3]

                        emitter['rotation'] = [w, x, y, z]

                        if fu.emitter_type == emitter_types[1][0]: # square emitter

                            emitter['type'] = 0
                            emitter['width'] = fu.object.dimensions[0]
                            emitter['height'] = fu.object.dimensions[1]

                        elif fu.emitter_type == emitter_types[2][0]: # circle emitter

                            emitter['type'] = 1
                            emitter['width'] = fu.object.dimensions[0] * 0.5
                            emitter['height'] = emitter['width']

                        emitters.append(emitter)
                        
            if len(fluidBlocks) > 0 or len(emitters) > 0 or '''len(geometry) > 0''':

                fluidModel = {
                    "viscosity": fm.viscosity,
                    "boundaryViscosity": fm.bviscosity,
                    "surfaceTension": fm.surften,
                    "adhesion": fm.adhesion,
                    "density0": fm.density0,
                }

                if len(fluidBlocks) > 0 :
                    fluidModel['fluidBlocks'] = fluidBlocks
                #if len(geometry) > 0 :
                    #fluidModel['geometry'] = []
                if len(emitters) > 0:
                    fluidModel['emitters'] = emitters
                    
                fluidSet['FluidModel'].append(fluidModel)

    if 'FluidModel' in fluidSet.keys():
        
        config['Fluid'] = {
        "viscosityMethod": viscosityMethod,
        "surfaceTensionMethod": surftenMethod,
        "adhesionMethod": adhesionMethod,
        "FluidModel": fluidSet['FluidModel']
    }

def export_boundary_set(object, config):

    boundarySetData = object.sph_boundary_set
    boundaryMethod = object.sph_boundary_set.boundaryMethod

    boundarySet = []

    for bm in boundarySetData.bs_list:

        if len(bm.bm_list) > 0:

            boundaryModel = {}

            if boundaryMethod == boundary_options[2][0]: # pcisph boundary handling
                boundaryModel['normalFct'] = bm.normalFct
                boundaryModel['tangFct'] = bm.tangFct

            boxes = []
            spheres = []
            geometries = []

            for bu in bm.bm_list:

                if bu.bu_type == boundary_unit_types[0][0]: # box

                    sideX = bu.object.dimensions[0]
                    sideY = bu.object.dimensions[1]
                    sideZ = bu.object.dimensions[2]

                    x = bu.object.location[0]
                    y = bu.object.location[1]
                    z = bu.object.location[2]

                    min = [x - sideX * 0.5,
                           y - sideY * 0.5,
                           z - sideZ * 0.5]
                    max = [x + sideX * 0.5,
                           y + sideY * 0.5,
                           z + sideZ * 0.5]

                    box = {
                        "min":min,
                        "max":max
                    }

                    if boundaryMethod == boundary_options[2][0]: # pcisph boundary handling
                        box['inverted'] = bu.inverted

                    boxes.append(box)

                elif bu.bu_type == boundary_unit_types[1][0]: # sphere

                    x = bu.object.location[0]
                    y = bu.object.location[1]
                    z = bu.object.location[2]

                    sphere = {
                        "pos": [x, y, z],
                        "radius": bu.object.dimensions[0] * 0.5
                    }

                    if boundaryMethod == boundary_options[2][0]: # pcisph boundary handling
                        sphere['inverted'] = bu.inverted

                    spheres.append(sphere)

                elif bu.bu_type == boundary_unit_types[2][0]: # geometry
                    
                    if bu.object is not None:
                        
                        bpy.ops.object.select_all(action='DESELECT')
                        bu.object.select_set(True)
                        
                        blend_file_path = bpy.data.filepath
                        directory = os.path.dirname(blend_file_path)
                        directory = os.path.join(directory, 'Models')
                        if os.path.isdir(directory) is False:
                            os.mkdir(directory)
                        target_file = os.path.join(directory, bu.object.name + '.obj')
                        
                        # Exportar obj
                        bpy.ops.export_scene.obj(filepath=target_file, use_selection=True, use_materials=False)
                        
                        geometry = {
                            "path": target_file,
                            "spacing": bu.spacing
                        }

                        geometries.append(geometry)

            if len(boxes) > 0 or len(spheres) > 0 or len(geometries) > 0:

                if len(boxes) > 0 :
                    boundaryModel['box'] = boxes
                if len(spheres) > 0:
                    boundaryModel['sphere'] = spheres
                if len(geometries) > 0 :
                    boundaryModel['geometry'] = geometries

                boundarySet.append(boundaryModel)


    if len(boundarySet) > 0:
        config['Boundary'] = []
        config['Boundary'] = boundarySet


###################################

def showMessageBox(title = "Message Box", icon = 'INFO', lines=""):
    def draw(self, context):
        for line in lines:
            self.layout.label(text=line)

    bpy.context.window_manager.popup_menu(draw, title = title, icon = icon)

# Funcion de update para los elementos de la interfaz
def set_modified(self,context):
    context.object.sph_modified = True
    
def updateMultiphase(self, context):
    fluidSet = context.object.sph_fluid_set
    
    # Dejar solo un fluid model cuando se desactiva multiphase 
    # (o cambiar el indice a 0 simplemente y hacer que solo se exporte el fm 0)
    if fluidSet.multiphase == False:
        if len(fluidSet.fs_list) > 1:
            for i in range(1, len(fluidSet.fs_list)):
                fluidSet.fs_list.remove(1)
        elif len(fluidSet.fs_list) == 0:
            fluidSet.fs_list.add()
            
        fluidSet.fm_index = 0
    
def deleteObject(name):
    
    # Deselect all
    bpy.ops.object.select_all(action='DESELECT')

    # Buscar el objeto para ver si existe
    obj = bpy.context.view_layer.objects.get(name)

    if obj is not None: # Si el objeto existe, eliminar su empty y luego el objeto
        if obj.hide_get() == True:
            obj.hide_set(False)
        obj.select_set(True)
        bpy.ops.object.delete()

def set_selected_fu(self, context):
    
    fluidSet = context.object.sph_fluid_set
    fm_index = fluidSet.fm_index
    fm = fluidSet.fs_list[fm_index]
    fu_index = fm.fu_index
    fm_list = fm.fm_list
    fu = fm_list[fu_index]
    
    if fu.object is not None:
        bpy.ops.object.select_all(action='DESELECT')
        fu.object.select_set(True)

def set_selected_bu(self, context):
    
    boundarySet = context.object.sph_boundary_set
    bm_index = boundarySet.bm_index
    bm = boundarySet.bs_list[bm_index]
    bu_index = bm.bu_index
    bm_list = bm.bm_list
    bu = bm_list[bu_index]
    
    if bu.object is not None:
        bpy.ops.object.select_all(action='DESELECT')
        bu.object.select_set(True)
        
def checkMeshBu(self, context):
    
    boundarySet = context.object.sph_boundary_set
    bm_index = boundarySet.bm_index
    bm = boundarySet.bs_list[bm_index]
    bu_index = bm.bu_index
    bm_list = bm.bm_list
    bu = bm_list[bu_index]
    
    print("checkmeshbu")
    if bu.object is not None:
        if bu.object.type != 'MESH':
            print("bu no mesh")
            showMessageBox("Warning", 'ERROR', ["Fluid unit must be a mesh"]) 
            bu.object = None
        
def checkMeshFu(self, context):
    
    fluidSet = context.object.sph_fluid_set
    fm_index = fluidSet.fm_index
    fm = fluidSet.fs_list[fm_index]
    fu_index = fm.fu_index
    fm_list = fm.fm_list
    fu = fm_list[fu_index]
    
    if fu.object is not None:
        if fu.object.type != 'MESH':
            showMessageBox("Warning", 'ERROR', ["Fluid unit must be a mesh"]) 
            fu.object = None
    

######### PROPERTIES ###################

class SceneInfo(bpy.types.PropertyGroup):

    startTime : bpy.props.FloatProperty(name="Start time",default=0.0, update=set_modified)
    endTime : bpy.props.FloatProperty(name="End time", default=5.0, update=set_modified)
    timeStep : bpy.props.FloatProperty(name="Time step", default=0.001, update=set_modified)
    fps : bpy.props.FloatProperty(name="FPS", default=120.0, update=set_modified)
    maxTimeStep : bpy.props.FloatProperty(name="Max. time step", default=0.05, update=set_modified)
    minTimeStep : bpy.props.FloatProperty(name="Min. time step", default=0.0, update=set_modified)
    particleRadius : bpy.props.FloatProperty(name="Particle radius", default=0.0, update=set_modified)
    gravity : bpy.props.FloatVectorProperty(name="Gravity", default=(0.0,0.0,-9.8), update=set_modified)

    simulationMethod : bpy.props.EnumProperty(
        name="Method",
        items=mode_options,
        description="Pressure solver",
        default="wcsph_method",
        update=set_modified
    )
    
    cflFactor : bpy.props.FloatProperty(name="CFL Factor", default=1.00, update=set_modified)

    # WCSPH
    stiffness : bpy.props.FloatProperty(name="Stiffness", default=100.0, update=set_modified)
    gamma : bpy.props.FloatProperty(name="Gamma", default=1.0, update=set_modified)
    
    # PCISPH
    minIterationsPCI : bpy.props.IntProperty(name="Min. iterations", default=3, update=set_modified)
    maxIterationsPCI : bpy.props.IntProperty(name="Max. iterations", default=100, update=set_modified)
    etaPCI : bpy.props.FloatProperty(name="Density error", default=0.01, update=set_modified)
    
    # DFSPH density
    minIterationsDF : bpy.props.IntProperty(name="Min. iterations", default=2, update=set_modified)
    maxIterationsDF : bpy.props.IntProperty(name="Max. iterations", default=100, update=set_modified)
    eta : bpy.props.FloatProperty(name="Density error", default=0.0005, update=set_modified)
    
    # DFSPH divergence
    divSolver : bpy.props.BoolProperty(name="Diverence solver", default=True, update=set_modified)
    minIterationsDFV : bpy.props.IntProperty(name="Min. iterations", default=1, update=set_modified)
    maxIterationsDFV : bpy.props.IntProperty(name="Max. iterations", default=100, update=set_modified)
    etaV : bpy.props.FloatProperty(name="Divergence error", default=0.001, update=set_modified)
    

class FluidUnit(PropertyGroup):
    """ Representa una unidad de fluido: fluid box, emitter o geometry """
    
    def clear(self):
        if self.object is not None and self.fu_type != fluid_unit_types[1][0]:
            deleteObject(self.object.name)
        if self.empty is not None:
            deleteObject(self.empty.name)
    
    object : bpy.props.PointerProperty(name="", type=bpy.types.Object, update=checkMeshFu)

    fu_type : bpy.props.EnumProperty(
        name="Fluid type",
        items=fluid_unit_types,
        description="Fluid unit type",
        default="fu_fluid_block_type",
        update=set_modified
    )

    # Todo lo que viene en la clase a partir de aqui es exclusivamente de los emisores, encontrar la manera de separarlo
    # Cuidado con crear primitivas al actualizar un enum (tambian se actualiza cuando se cambia el valor desde el codigo ademas de en la interfaz)
    def update_emitter_enum(self, context):

        if self.emitter_type != emitter_types[0][0]:
            
            # guardar objeto activo antes de crear el emisor y el empty
            nameActiveObject = bpy.context.view_layer.objects.active.name
            
            # Añadir un empty para indicar la dirección del emisor
            bpy.ops.object.empty_add(type='SINGLE_ARROW', align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
            self.empty = bpy.context.view_layer.objects.active
            
            if self.emitter_type == emitter_types[1][0]: # Square

                # crear el plano
                bpy.ops.mesh.primitive_plane_add(size=1, enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(0.5, 0.5, 0.5))
                self.object = bpy.context.view_layer.objects.active
                self.object.name = "Square emitter"

            if self.emitter_type == emitter_types[2][0]: # Circle

                # crear el circulo y almacenarlo
                bpy.ops.mesh.primitive_circle_add(radius=0.5, enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(0.5, 0.5, 0.5))
                self.object = bpy.context.view_layer.objects.active
                self.object.name = "Circle emitter"

            # hacer al empty hijo del objeto emisor
            self.empty.parent = self.object

            # reestablecer el objeto activo
            bpy.context.view_layer.objects.active = bpy.data.objects[nameActiveObject]

    emitter_type : bpy.props.EnumProperty(
        name="Emitter type",
        items=emitter_types,
        description="Fluid unit type",
        default="none_emitter",
        update=update_emitter_enum
    )
    
    empty : bpy.props.PointerProperty(type=bpy.types.Object)

    velocity : bpy.props.FloatProperty(name="Velocity", default=0.2, update=set_modified)
    spacing : bpy.props.FloatProperty(name="Spacing", default=1.0, update=set_modified)
    numParticles: bpy.props.IntProperty(name="Number of particles", default=10000, update=set_modified)
    startTime : bpy.props.FloatProperty(name="Emit start time", default=0.0, update=set_modified)


class FluidModel(PropertyGroup):
    """Group of properties representing fluid model in the fluid Set list."""

    def clear(self):
        for fu in self.fm_list:
            fu.clear()            

    name : bpy.props.StringProperty(name="", default="Fluid model", update=set_modified) # Cuidado no tiene nombre
    density0 : bpy.props.FloatProperty(name="Reference density", default=1000.0, update=set_modified)
    viscosity : bpy.props.FloatProperty(name="Viscosity", default=0.01, update=set_modified)
    bviscosity : bpy.props.FloatProperty(name="Boundary viscosity", default=0.0, update=set_modified)
    surften : bpy.props.FloatProperty(name="Surface Tension", default=0.15, update=set_modified)
    adhesion : bpy.props.FloatProperty(name="Adhesion", default=1.0, update=set_modified)

    fm_list : CollectionProperty(type = FluidUnit)
    fu_index : bpy.props.IntProperty(name="Fluid unit index", default=0, update=set_selected_fu)

    current_fu_type : bpy.props.EnumProperty(
        name="Fluid type",
        items=fluid_unit_types,
        description="Fluid unit type",
        default="fu_fluid_block_type",
        update=set_modified
    )

class BoundaryUnit(PropertyGroup):
    
    def clear(self):
        if self.object is not None and self.bu_type != boundary_unit_types[2][0]:
            deleteObject(self.object.name)

    # necesarias en pcisph boundary handling
    inverted: bpy.props.BoolProperty(name="Inverted", default=False, update=set_modified)

    # necesario para la geometria
    spacing : bpy.props.FloatProperty(name="Spacing", default=0.7, update=set_modified)

    object : bpy.props.PointerProperty(name="", type=bpy.types.Object, update=checkMeshBu)

    bu_type : bpy.props.EnumProperty(
        name="Boundary type",
        items=boundary_unit_types,
        description="Boundary unit type",
        default="bu_box_type",
        update=set_modified
    )



class BoundaryModel(PropertyGroup):

    def clear(self):
        for bu in self.bm_list:
            bu.clear()

    bm_list : CollectionProperty(type = BoundaryUnit)
    bu_index : bpy.props.IntProperty(name="Boundary unit index", default=0, update=set_selected_bu)
    name : bpy.props.StringProperty(name="", default="Boundary model", update=set_modified) # Cuidado no tiene nombre

    # necesario en pcisph boundary handling
    normalFct : bpy.props.FloatProperty(name="Normal conservation", default=0.0, update=set_modified)
    tangFct : bpy.props.FloatProperty(name="Tangencial conservation", default=1.0, update=set_modified)

    current_bu_type : bpy.props.EnumProperty(
        name="Boundary type",
        items=boundary_unit_types,
        description="Boundary unit type",
        default="bu_box_type",
        update=set_modified
    )



class BoundarySet(bpy.types.PropertyGroup):
    
    def clear(self):
        for bm in self.bs_list:
            bm.clear()

    bs_list : CollectionProperty(type = BoundaryModel)
    bm_index : bpy.props.IntProperty(name="Boundary model index", default=0, update=set_modified)

    boundaryMethod : bpy.props.EnumProperty(
        name="Boundary method",
        items=boundary_options,
        description="Boundary handling Method",
        default="akinci_boundary_method",
        update=set_modified
    )


class FluidSet(bpy.types.PropertyGroup):
    
    def clear(self):
        for fm in self.fs_list:
            fm.clear()

    viscosityMethod : bpy.props.EnumProperty(
        name="Viscosity Method",
        items=viscosity_options,
        description="Viscosity Method",
        default="artificial_visco_method",
        update=set_modified
    )

    surftenMethod : bpy.props.EnumProperty(
        name="Surface Tension Method",
        items=surften_options,
        description="Surface Tension Method",
        default="akinci_st_method",
        update=set_modified
    )

    adhesionMethod : bpy.props.EnumProperty(
        name="Adhesion Method",
        items=adhesion_options,
        description="Adhesion Method",
        default="none_adhesion_method",
        update=set_modified
    )

    fs_list : CollectionProperty(name="List of fluid models", type = FluidModel)
    fm_index : bpy.props.IntProperty(name="Fluid model index", default=0, update=set_modified)
    multiphase : bpy.props.BoolProperty(name="Multiphase simulation", default=False, update=updateMultiphase)

##### Fluid Set UIList

class FluidSetList(UIList):
    bl_idname = "OBJECT_UL_fluidsetlist"

    def draw_item(self, context, layout, data, item, icon, active_data,
                  active_propname, index):

        custom_icon = 'MATFLUID'

        # Make sure your code supports all 3 layout types
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            split = layout.split()
            split.prop(item, "name", emboss=False, icon=custom_icon)

        elif self.layout_type in {'GRID'}:
            layout.alignment = 'CENTER'
            split = layout.split()
            split.prop(item, "name", emboss=False, icon=custom_icon)

#### Fluid Model UIList

class FluidModelList(UIList):
    bl_idname = "OBJECT_UL_fluidmodellist"

    def draw_item(self, context, layout, data, item, icon, active_data,
                  active_propname, index):
        custom_icon = 'UNLINKED'
        # comprobar si sirve de algo para cuando el objeto es none
        if item.object is not None:
            if item.fu_type == fluid_unit_types[0][0]: # fluid block
                custom_icon = 'CUBE'
            
            elif item.fu_type == fluid_unit_types[1][0]: # Geometry
                custom_icon = 'MESH_MONKEY'

            elif item.fu_type == fluid_unit_types[2][0]: # Emitter
                
                if item.emitter_type == emitter_types[1][0]: #square
                    custom_icon = 'MESH_PLANE'
                elif item.emitter_type == emitter_types[2][0]: #circle
                    custom_icon = 'MESH_CIRCLE'
                    
            split = layout.split()
            split.prop(item.object, "name", text="", emboss=False, icon=custom_icon)

        else:

            split = layout.split()

            if item.fu_type == fluid_unit_types[1][0]: # Geometry
                split.label(text="Empty geometry", icon=custom_icon)

            elif item.fu_type == fluid_unit_types[2][0]: # Emitter                
                split.label(text="Empty emitter", icon=custom_icon)

##### Boundary Set UIList
class BoundarySetList(UIList):
    bl_idname = "OBJECT_UL_boundarysetlist"

    def draw_item(self, context, layout, data, item, icon, active_data,
                  active_propname, index):

        # We could write some code to decide which icon to use here...
        custom_icon = 'RIGID_BODY'

        # Make sure your code supports all 3 layout types
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            split = layout.split()
            split.prop(item, "name", emboss=False, icon=custom_icon)

        elif self.layout_type in {'GRID'}:
            layout.alignment = 'CENTER'
            split = layout.split()
            split.prop(item, "name", emboss=False, icon=custom_icon)

##### Boundary Model UIList
class BoundaryModelList(UIList):
    bl_idname = "OBJECT_UL_boundarymodellist"

    def draw_item(self, context, layout, data, item, icon, active_data,
                  active_propname, index):
        custom_icon = 'UNLINKED'

        if item.object is not None:

            if item.bu_type == boundary_unit_types[0][0]: # box
                custom_icon = 'CUBE'
            elif item.bu_type == boundary_unit_types[1][0]: # Sphere
                custom_icon = 'SPHERE'
            elif item.bu_type == boundary_unit_types[2][0]: # geometry
                custom_icon = 'MESH_MONKEY'

            split = layout.split()
            split.prop(item.object, "name", text="", emboss=False, icon=custom_icon)

        else:

            split = layout.split()

            if item.bu_type == boundary_unit_types[2][0]: # geometry
                split.label(text="Empty geometry", icon=custom_icon)

####### PANELES ##############

# Panel del que heredan todos los demas
class MyPanel(bpy.types.Panel):
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "physics"

# Panel principal que contiene los demas panels
class MainPanel(MyPanel):
    bl_label = "SPH Fluid"
    bl_idname = "OBJECT_PT_sph_main"
    bl_icon = "ADD"

    def draw(self, context):
        obj = context.object
        layout = self.layout

        if obj.sph_created is False:
            split = layout.split()

            col = split.column()
            col.operator("sph.createscene", icon="ADD")

            col = split.column()
            col.operator("sph.importscene", icon="IMPORT")
        else:
            split = layout.split()
            col = split.column();
            col.operator("sph.deletescene", icon="TRASH")

            col = split.column()
            col.operator("sph.exportscene", icon="EXPORT")
            
            row = layout.row()
            row.operator("sph.readframes", icon="RENDER_ANIMATION")

# Subpanel para la configuracion general
class GeneralConfPanel(MyPanel):
    bl_label = "General configuration"
    bl_parent_id = "OBJECT_PT_sph_main"
    bl_idname = "OBJECT_PT_sph_general_conf"

    @classmethod
    def poll(self, context):
        return context.object.sph_created

    def draw(self, context):
            obj = context.object
            layout = self.layout
            #layout.use_property_split = True # Otra distribucion
            #layout.use_property_decorate = False

            if obj.sph_created is True:

                split = layout.split()

                col = split.column()
                col.prop(obj.sph_scene_data, "startTime")

                col = split.column()
                col.prop(obj.sph_scene_data, "endTime")

                split = layout.split()

                col = split.column()
                col.prop(obj.sph_scene_data, "timeStep")

                col = split.column()
                col.prop(obj.sph_scene_data, "fps")

                split = layout.split()

                col = split.column()
                col.prop(obj.sph_scene_data, "minTimeStep")

                col = split.column()
                col.prop(obj.sph_scene_data, "maxTimeStep")

                row = layout.row()
                row.prop(obj.sph_scene_data, "particleRadius")

                row = layout.row()
                row.prop(obj.sph_scene_data, "gravity")

# Subpanel para los parametros del solver
class solverSubPanel(MyPanel):
    bl_label = "Solver"
    bl_parent_id = "OBJECT_PT_sph_main"
    bl_idname = "OBJECT_PT_sph_solver"

    @classmethod
    def poll(self, context):
        return context.object.sph_created

    def draw(self, context):
        obj = context.object
        layout = self.layout
        layout.use_property_decorate = False

        if obj.sph_created is True:
            row = layout.row()
            row.prop(obj.sph_scene_data, "simulationMethod")
            simulationMethod = obj.sph_scene_data.simulationMethod
            layout.use_property_split = True # Otra distribucion
            if simulationMethod == "wcsph_method":

                row = layout.row()
                row.prop(obj.sph_scene_data, "stiffness")

                row = layout.row()
                row.prop(obj.sph_scene_data, "gamma")

                row = layout.row()
                row.prop(obj.sph_scene_data, "cflFactor")

            elif  simulationMethod == "pcisph_method":

                row = layout.row()
                row.prop(obj.sph_scene_data, "etaPCI")
                
                row = layout.row()
                row.prop(obj.sph_scene_data, "minIterationsPCI")
                
                row = layout.row()
                row.prop(obj.sph_scene_data, "maxIterationsPCI")

                row = layout.row()
                row.prop(obj.sph_scene_data, "cflFactor")

            elif  simulationMethod == "dfsph_method":

                row = layout.row()
                row.prop(obj.sph_scene_data, "eta")
                
                row = layout.row()
                row.prop(obj.sph_scene_data, "minIterationsDF")
                
                row = layout.row()
                row.prop(obj.sph_scene_data, "maxIterationsDF")

                row = layout.row()
                row.prop(obj.sph_scene_data, "cflFactor")
                
                row = layout.row()
                row.prop(obj.sph_scene_data, "divSolver")
                
                if obj.sph_scene_data.divSolver == True:
                    
                    row = layout.row()
                    row.prop(obj.sph_scene_data, "etaV")
                    
                    row = layout.row()
                    row.prop(obj.sph_scene_data, "minIterationsDFV")
                    
                    row = layout.row()
                    row.prop(obj.sph_scene_data, "maxIterationsDFV")
                    
                    

# Sub panel para el fluid set
class fluidSetSubPanel(MyPanel):
    bl_label = "Fluid Set"
    bl_parent_id = "OBJECT_PT_sph_main"
    bl_idname = "OBJECT_PT_sph_fluid_set"

    @classmethod
    def poll(self, context):
        return context.object.sph_created

    def draw(self, context):
        obj = context.object
        layout = self.layout
        layout.use_property_split = True # Otra distribucion
        layout.use_property_decorate = False

        if obj.sph_created is True:
            row = layout.row()
            row.prop(obj.sph_fluid_set, "viscosityMethod")
            row = layout.row()
            row.prop(obj.sph_fluid_set, "surftenMethod")
            row = layout.row()
            row.prop(obj.sph_fluid_set, "adhesionMethod")
            row = layout.row()
            row.prop(obj.sph_fluid_set, "multiphase")
            
            if obj.sph_fluid_set.multiphase:
                
                # FluidSet of Fluid Models
                row = layout.row()
                row.template_list("OBJECT_UL_fluidsetlist", "The_List", obj.sph_fluid_set,
                                  "fs_list", obj.sph_fluid_set, "fm_index")
       
                col = row.column(align=True)
                col.operator("sph.addfluidmodel", icon="ADD")
                col.operator("sph.removefluidmodel", icon="REMOVE")


            if obj.sph_fluid_set.fm_index >= 0 and obj.sph_fluid_set.fs_list:
                fm = obj.sph_fluid_set.fs_list[obj.sph_fluid_set.fm_index]
                
                row = layout.row()
                row.prop(fm, "density0")

                row = layout.row()
                row.prop(fm, "viscosity")

                row = layout.row()
                row.prop(fm, "bviscosity")

                row = layout.row()
                row.prop(fm, "surften")

                row = layout.row()
                row.prop(fm, "adhesion")

                # FluidModel of Fluid Units                                  
                row = layout.row()
                col = row.column()
                col.template_list("OBJECT_UL_fluidmodellist", "The_List", fm,
                              "fm_list", fm, "fu_index")
                              
                col = col.row()
                col.use_property_split = False
                col.prop(fm, "current_fu_type", expand=True)
                              
                col = row.column(align=True)
                col.operator("sph.addfluidunit", icon="ADD")
                col.operator("sph.removefluidunit", icon="REMOVE")

                if fm.fu_index >= 0 and fm.fm_list:
                    fu = fm.fm_list[fm.fu_index]

                    if fu.fu_type == fluid_unit_types[0][0]: # fluid block

                        row = layout.row()
                        row.prop(fu.object, "location")

                        row = layout.row()
                        row.prop(fu.object, "rotation_euler")

                        row = layout.row()
                        row.prop(fu.object, "scale")

                    elif fu.fu_type == fluid_unit_types[1][0]: # geometry

                        layout.use_property_split = False
                        row = layout.row()
                        row.prop_search(fu, "object", context.scene, "objects")
                        layout.use_property_split = True

                        if fu.object is not None:

                            row = layout.row()
                            row.prop(fu.object, "location")

                            row = layout.row()
                            row.prop(fu.object, "rotation_euler")

                            row = layout.row()
                            row.prop(fu.object, "scale")

                    elif fu.fu_type == fluid_unit_types[2][0]: # emitter

                        if fu.object is None:
                            row = layout.row()
                            row.prop(fu, "emitter_type")

                            row = layout.row()
                            row.prop(fu, "startTime")

                            row = layout.row()
                            row.prop(fu, "velocity")

                            row = layout.row()
                            row.prop(fu, "numParticles")
                            
                            row = layout.row()
                            row.prop(fu, "spacing")

                        else:

                            row = layout.row()
                            row.prop(fu, "startTime")

                            row = layout.row()
                            row.prop(fu, "velocity")

                            row = layout.row()
                            row.prop(fu, "numParticles")
                            
                            row = layout.row()
                            row.prop(fu, "spacing")

                            row = layout.row()
                            row.prop(fu.object, "location")

                            row = layout.row()
                            row.prop(fu.object, "rotation_euler")

                            row = layout.row()
                            row.prop(fu.object, "scale")

# Subpanel para el Boundary set
class boundarySetSubPanel(MyPanel):
    bl_label = "Boundary Set"
    bl_parent_id = "OBJECT_PT_sph_main"
    bl_idname = "OBJECT_PT_sph_boundary_set"

    @classmethod
    def poll(self, context):
        return context.object.sph_created

    def draw(self, context):
        obj = context.object
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False
        

        if obj.sph_created is True:
            row = layout.row()
            row.prop(obj.sph_boundary_set, "boundaryMethod")

            if obj.sph_boundary_set.boundaryMethod != boundary_options[0][0]:

                row = layout.row()
                row.template_list("OBJECT_UL_boundarysetlist", "The_List", obj.sph_boundary_set,
                              "bs_list", obj.sph_boundary_set, "bm_index")

                col = row.column(align=True)
                col.operator("sph.addboundarymodel", icon="ADD")
                col.operator("sph.removeboundarymodel", icon="REMOVE")

            if obj.sph_boundary_set.bm_index >= 0 and obj.sph_boundary_set.bs_list:
                bm = obj.sph_boundary_set.bs_list[obj.sph_boundary_set.bm_index]

                row = layout.row()
                col = row.column()
                col.template_list("OBJECT_UL_boundarymodellist", "The_List", bm,
                              "bm_list", bm, "bu_index")
                              
                col = col.row()
                col.use_property_split = False
                col.prop(bm, "current_bu_type", expand=True)
                              
                col = row.column(align=True)
                col.operator("sph.addboundaryunit", icon="ADD")
                col.operator("sph.removeboundaryunit", icon="REMOVE")

                if bm.bu_index >= 0 and bm.bm_list:
                    bu = bm.bm_list[bm.bu_index]

                    if bu.bu_type == boundary_unit_types[0][0]: # box

                        if obj.sph_boundary_set.boundaryMethod == boundary_options[2][0]:
                            row = layout.row()
                            row.prop(bu, "inverted")

                        row = layout.row()
                        row.prop(bu.object, "location")

                        row = layout.row()
                        row.prop(bu.object, "rotation_euler")

                        row = layout.row()
                        row.prop(bu.object, "scale")

                    elif bu.bu_type == boundary_unit_types[1][0]: # sphere
                        
                        if obj.sph_boundary_set.boundaryMethod == boundary_options[2][0]:
                            row = layout.row()
                            row.prop(bu, "inverted")

                        row = layout.row()
                        row.prop(bu.object, "location")

                        row = layout.row()
                        row.prop(bu.object, "rotation_euler")

                        row = layout.row()
                        row.prop(bu.object, "scale")

                    elif bu.bu_type == boundary_unit_types[2][0]: # geometry

                        layout.use_property_split = False
                        row = layout.row()
                        row.prop_search(bu, "object", context.scene, "objects")
                        layout.use_property_split = True

                        if bu.object is not None:

                            row = layout.row()
                            row.prop(bu, "spacing")

                            row = layout.row()
                            row.prop(bu.object, "location")

                            row = layout.row()
                            row.prop(bu.object, "rotation_euler")

                            row = layout.row()
                            row.prop(bu.object, "scale")


############################### OPERADORES
class readFramesOp(bpy.types.Operator, ImportHelper):
    bl_idname = "sph.readframes"
    bl_label = "Read frames"
    
    files : CollectionProperty(name='File paths', type=bpy.types.OperatorFileListElement)
    
    def execute(self, context):
        read_frames(context, self.files, self.filepath)
        
        return {"FINISHED"}

class createSceneOp(bpy.types.Operator):
    bl_idname = "sph.createscene"
    bl_label = "Create new scene"

    def execute(self, context):
        obj = context.object
        obj.sph_created = True
        
        obj.sph_fluid_set.fs_list.add()
        obj.sph_fluid_set.fm_index = 0
        obj.display_type = 'WIRE'

        return {"FINISHED"}

class importSceneOp(bpy.types.Operator, ImportHelper):
    bl_idname = "sph.importscene"
    bl_label = "Import scene"

    def execute(self, context):
        obj = context.object
        obj.sph_created = True
        obj.display_type = 'WIRE'
        
        import_config(obj, self.filepath)

        return {"FINISHED"}

class exportSceneOp(bpy.types.Operator):
    bl_idname = "sph.exportscene"
    bl_label = "Export scene"

    def execute(self, context):

        return export_config(context.object)

class deleteSceneOp(bpy.types.Operator):
    bl_idname = "sph.deletescene"
    bl_label = "Delete scene"

    def execute(self, context):
        obj = context.object
        obj.sph_created = False

        # Resetear fluid set y boundary set
        # Eliminar objetos asociados
        obj.sph_fluid_set.clear() 
        obj.sph_boundary_set.clear() 
        obj.sph_fluid_set.multiphase = False
        
        # no borra los objetos asociados a los fluid units
        obj.sph_fluid_set.fs_list.clear()    
        obj.sph_boundary_set.bs_list.clear() 
        
        obj.display_type = 'TEXTURED'

        return {"FINISHED"}

class addFluidModelOp(bpy.types.Operator):
    bl_idname = "sph.addfluidmodel"
    bl_label = ""

    def execute(self, context):
        fs_list = context.object.sph_fluid_set.fs_list
        fm_index = context.object.sph_fluid_set.fm_index
        
        new_index = len(fs_list)
        fs_list.add()
        context.object.sph_fluid_set.fm_index = new_index

        return{'FINISHED'}


# Borrar fluidModel
class removeFluidModelOp(Operator):

    bl_idname = "sph.removefluidmodel"
    bl_label = ""

    @classmethod
    def poll(cls, context):
        return context.object.sph_fluid_set.fs_list

    def execute(self, context):
        fs_list = context.object.sph_fluid_set.fs_list
        index = context.object.sph_fluid_set.fm_index
        fm = fs_list[index]

        # Eliminar objetos asociados
        fm.clear()
                
        fs_list.remove(index)
        context.object.sph_fluid_set.fm_index = min(max(0, index - 1), len(fs_list) - 1)

        return{'FINISHED'}


# Añadir fluid unit
class addFluidUnitOp(bpy.types.Operator):
    bl_idname = "sph.addfluidunit"
    bl_label = ""

    def execute(self, context):
        fm_index = context.object.sph_fluid_set.fm_index
        fm = context.object.sph_fluid_set.fs_list[fm_index]
        fm.fm_list.add()
        fm.fu_index = len(fm.fm_list) - 1

        if fm.current_fu_type == fluid_unit_types[0][0]: # fluid block
            # guardar objeto activo antes de crear el nuevo cubo
            nameActiveObject = bpy.context.view_layer.objects.active.name

            # crear el cubo
            bpy.ops.mesh.primitive_cube_add(size=0.5, enter_editmode=False, align="WORLD", location=(0, 0, 0), scale=(1, 1, 1))
            bpy.context.view_layer.objects.active.name = "Fluid block"

            # Guardar en el fluid unit una referencia al objeto que lo representa y su tipo
            fm.fm_list[len(fm.fm_list) - 1].object = bpy.context.view_layer.objects.active
            fm.fm_list[len(fm.fm_list) - 1].fu_type = fm.current_fu_type

            # reestablecer el objeto activo
            bpy.context.view_layer.objects.active = bpy.data.objects[nameActiveObject]


        elif fm.current_fu_type == fluid_unit_types[1][0]: # geometry
            fm.fm_list[len(fm.fm_list) - 1].fu_type = fm.current_fu_type
        elif fm.current_fu_type == fluid_unit_types[2][0]: # emmiter
            fm.fm_list[len(fm.fm_list) - 1].fu_type = fm.current_fu_type

        return{'FINISHED'}


# Borrar fluidUnit
class removeFluidUnitOp(Operator):

    bl_idname = "sph.removefluidunit"
    bl_label = ""

    @classmethod
    def poll(cls, context):
        fm = context.object.sph_fluid_set.fs_list[context.object.sph_fluid_set.fm_index]
        return fm.fm_list

    def execute(self, context):
        fm = context.object.sph_fluid_set.fs_list[context.object.sph_fluid_set.fm_index]
        fm_list = fm.fm_list
        index = fm.fu_index

        # Eliminar objetos asociados
        fm_list[index].clear()

        # borrar elemento de la lista
        fm_list.remove(index)
        fm.fu_index = min(max(0, index - 1), len(fm_list) - 1)

        return{'FINISHED'}

# Añadir boundary unit
class addBoundaryUnitOp(bpy.types.Operator):
    bl_idname = "sph.addboundaryunit"
    bl_label = ""

    def execute(self, context):
        bm_index = context.object.sph_boundary_set.bm_index
        bm = context.object.sph_boundary_set.bs_list[bm_index]
        bm.bm_list.add()
        bm.bu_index = len(bm.bm_list) - 1

        if bm.current_bu_type == boundary_unit_types[0][0]: # box
             # guardar objeto activo antes de crear el nuevo cubo
            nameActiveObject = bpy.context.view_layer.objects.active.name

            # crear el cubo
            bpy.ops.mesh.primitive_cube_add(size=0.5, enter_editmode=False, align="WORLD", location=(0, 0, 0), scale=(1, 1, 1))
            bpy.context.view_layer.objects.active.name = "Boundary cube"

            # Guardar en el fluid unit una referencia al objeto que lo representa y su tipo
            bm.bm_list[len(bm.bm_list) - 1].object = bpy.context.view_layer.objects.active
            bm.bm_list[len(bm.bm_list) - 1].bu_type = bm.current_bu_type

            # reestablecer el objeto activo
            bpy.context.view_layer.objects.active = bpy.data.objects[nameActiveObject]

        elif bm.current_bu_type == boundary_unit_types[1][0]: # sphere
            # guardar objeto activo antes de crear el nuevo cubo
            nameActiveObject = bpy.context.view_layer.objects.active.name

            # crear el cubo
            bpy.ops.mesh.primitive_uv_sphere_add(radius=0.5, enter_editmode=False, align="WORLD", location=(0, 0, 0), scale=(1, 1, 1))
            bpy.context.view_layer.objects.active.name = "Boundary sphere"

            # Guardar en el fluid unit una referencia al objeto que lo representa y su tipo
            bm.bm_list[len(bm.bm_list) - 1].object = bpy.context.view_layer.objects.active
            bm.bm_list[len(bm.bm_list) - 1].bu_type = bm.current_bu_type

            # reestablecer el objeto activo
            bpy.context.view_layer.objects.active = bpy.data.objects[nameActiveObject]

        elif bm.current_bu_type == boundary_unit_types[2][0]: # geometry
            bm.bm_list[len(bm.bm_list) - 1].bu_type = bm.current_bu_type

        return{'FINISHED'}


# Borrar boundaryUnit
class removeBoundaryUnitOp(Operator):

    bl_idname = "sph.removeboundaryunit"
    bl_label = ""

    @classmethod
    def poll(cls, context):
        bm = context.object.sph_boundary_set.bs_list[context.object.sph_boundary_set.bm_index]
        return bm.bm_list

    def execute(self, context):
        bm = context.object.sph_boundary_set.bs_list[context.object.sph_boundary_set.bm_index]
        bm_list = bm.bm_list
        index = bm.bu_index

        # Eliminar objetos asociados
        bm_list[index].clear()

        # borrar elemento de la lista
        bm_list.remove(index)
        bm.bu_index = min(max(0, index - 1), len(bm_list) - 1)

        return{'FINISHED'}


# Añadir boundadary model
class addBoundaryModelOp(bpy.types.Operator):
    bl_idname = "sph.addboundarymodel"
    bl_label = ""

    def execute(self, context):
        bs_list = context.object.sph_boundary_set.bs_list
        bs_list.add()
        context.object.sph_boundary_set.bm_index = len(bs_list) - 1

        return{'FINISHED'}


# Borrar boundary Model
class removeBoundaryModelOp(Operator):

    bl_idname = "sph.removeboundarymodel"
    bl_label = ""

    @classmethod
    def poll(cls, context):
        return context.object.sph_boundary_set.bs_list

    def execute(self, context):
        bs_list = context.object.sph_boundary_set.bs_list
        index = context.object.sph_boundary_set.bm_index
        bm = bs_list[index]
        
        bm.clear()

        bs_list.remove(index)
        context.object.sph_boundary_set.bm_index = min(max(0, index - 1), len(bs_list) - 1)

        return{'FINISHED'}

def register():
    bpy.utils.register_class(MainPanel)
    bpy.utils.register_class(GeneralConfPanel)
    bpy.utils.register_class(solverSubPanel)
    bpy.utils.register_class(fluidSetSubPanel)
    bpy.utils.register_class(boundarySetSubPanel)

    bpy.utils.register_class(SceneInfo)
    bpy.utils.register_class(FluidUnit)
    bpy.utils.register_class(FluidModel)
    bpy.utils.register_class(FluidSet)
    bpy.utils.register_class(BoundaryUnit)
    bpy.utils.register_class(BoundaryModel)
    bpy.utils.register_class(BoundarySet)

    bpy.utils.register_class(createSceneOp)
    bpy.utils.register_class(importSceneOp)
    bpy.utils.register_class(exportSceneOp)
    bpy.utils.register_class(deleteSceneOp)
    bpy.utils.register_class(readFramesOp)
    bpy.utils.register_class(addFluidModelOp)
    bpy.utils.register_class(removeFluidModelOp)
    bpy.utils.register_class(addBoundaryModelOp)
    bpy.utils.register_class(removeBoundaryModelOp)
    bpy.utils.register_class(addFluidUnitOp)
    bpy.utils.register_class(removeFluidUnitOp)
    bpy.utils.register_class(addBoundaryUnitOp)
    bpy.utils.register_class(removeBoundaryUnitOp)

    bpy.utils.register_class(FluidSetList)
    bpy.utils.register_class(FluidModelList)
    bpy.utils.register_class(BoundarySetList)
    bpy.utils.register_class(BoundaryModelList)

    bpy.types.Object.sph_modified = bpy.props.BoolProperty(name = "SPH Scene modified", default = False)
    bpy.types.Object.sph_created = bpy.props.BoolProperty(name = "SPH Scene created", default = False)
    bpy.types.Object.sph_scene_data = bpy.props.PointerProperty(type=SceneInfo)
    bpy.types.Object.sph_fluid_set = bpy.props.PointerProperty(type=FluidSet)
    bpy.types.Object.sph_boundary_set = bpy.props.PointerProperty(type=BoundarySet)
    

def unregister():
    bpy.utils.unregister_class(MainPanel)
    bpy.utils.unregister_class(GeneralConfPanel)
    bpy.utils.unregister_class(solverSubPanel)
    bpy.utils.unregister_class(fluidSetSubPanel)
    bpy.utils.unregister_class(boundarySetSubPanel)

    bpy.utils.unregister_class(SceneInfo)
    bpy.utils.unregister_class(FluidUnit)
    bpy.utils.unregister_class(FluidModel)
    bpy.utils.unregister_class(FluidSet)
    bpy.utils.unregister_class(BoundaryUnit)
    bpy.utils.unregister_class(BoundaryModel)
    bpy.utils.unregister_class(BoundarySet)

    bpy.utils.unregister_class(createSceneOp)
    bpy.utils.unregister_class(importSceneOp)
    bpy.utils.unregister_class(exportSceneOp)
    bpy.utils.unregister_class(deleteSceneOp)
    bpy.utils.unregister_class(readFramesOp)
    bpy.utils.unregister_class(addFluidModelOp)
    bpy.utils.unregister_class(removeFluidModelOp)
    bpy.utils.unregister_class(addBoundaryModelOp)
    bpy.utils.unregister_class(removeBoundaryModelOp)
    bpy.utils.unregister_class(addFluidUnitOp)
    bpy.utils.unregister_class(removeFluidUnitOp)
    bpy.utils.unregister_class(addBoundaryUnitOp)
    bpy.utils.unregister_class(removeBoundaryUnitOp)

    bpy.utils.unregister_class(FluidSetList)
    bpy.utils.unregister_class(FluidModelList)
    bpy.utils.unregister_class(BoundarySetList)
    bpy.utils.unregister_class(BoundaryModelList)

    del bpy.types.Object.sph_scene_data
    del bpy.types.Object.sph_fluid_set
    del bpy.types.Object.sph_boundary_set
    del bpy.types.Object.sph_modified
    del bpy.types.Object.sph_created


if __name__ == "__main__" :
    register()



