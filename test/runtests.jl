using Test
using SimplSurfs


@testset "Quad-Edge Data Structure Tests" begin

    @testset "make_edge: Pointer Consistency" begin
        mesh = PolyhedralMesh{0,Any}()
        e = make_edge(mesh)

        # Test rotational pointers
        @test rot(rot(rot(rot(e)))) === e
        @test rot(e) !== e
        @test rot(rot(e)) !== e

        # Test flip pointers
        @test flip(flip(e)) === e
        @test flip(e) !== e
        @test rot(rot(flip(e))) === flip(rot(rot(e))) # flip and rot commute for 180 degrees

        # Test initial 'next' pointers for an isolated edge
        @test next(e) === e
        @test next(rot(rot(e))) === rot(rot(e))
        @test next(rot(e)) === rot(rot(rot(e)))
        @test next(rot(rot(rot(e)))) === rot(e)
    end

    @testset "Data Assignment and Accessors" begin
        mesh = PolyhedralMesh{0,Any}()
        e = make_edge(mesh)

        v1 = Vertex(1, mesh=mesh)
        v2 = Vertex(2, mesh=mesh)
        f1 = Face(1, mesh=mesh)
        f2 = Face(2, mesh=mesh)

        # Set data using the mutator functions
        tail!(e, v1)
        head!(e, v2)
        left!(e, f1)
        right!(e, f2)

        # Test accessors
        @test tail(e) === v1
        @test head(e) === v2
        @test left(e) === f1
        @test right(e) === f2

        # Test propagation to other edges in the group
        @test tail(flip(e)) === v1 # Flipped edge shares the same tail
        @test head(flip(e)) === v2 # and head
        @test left(flip(e)) === f2 # but left/right are swapped
        @test right(flip(e)) === f1

        @test tail(rot(e)) === f2 # right(e) is tail(rot(e))
        @test head(rot(e)) === f1 # left(e) is head(rot(e))
    end

    @testset "splice!: Joining and Splitting Rings" begin
        mesh = PolyhedralMesh{0,Any}()
        v1 = Vertex(1, mesh=mesh)
        v2 = Vertex(2, mesh=mesh)

        a = make_edge(mesh)
        b = make_edge(mesh)
        tail!(a, v1)
        tail!(b, v2)

        # --- Test Joining ---
        # Initially, a and b are in separate 1-element rings
        @test next(a) === a
        @test next(b) === b
        @test next(flip(a)) === flip(a)
        @test next(flip(b)) === flip(b)

        # Manually unify tails to test ring joining
        # tail!(b, v1)

        # Splice with unified tails.
        SimplSurfs.splice!(a, b)

        @test next(a) === b
        @test next(b) === a
        @test next(flip(a)) === flip(b)
        @test next(flip(b)) === flip(a)


        # --- Test Splitting ---
        # Now split the 2-element ring we just created
        SimplSurfs.splice!(a, b)
        @test next(a) === a
        @test next(b) === b
        @test next(flip(a)) === flip(a)
        @test next(flip(b)) === flip(b)
    end
end



@testset "Polyhedral Mesh Tests" begin

    @testset "Polyhedral Mesh Construction" begin

        """
            _check_consistency(mesh::PolyhedralMesh)

        Check the consistency of the mesh, i.e. for each edge `e`: 
        - tail(next(e)) === tail(e) 
        - left(edge) === right(next(e))
        - tail(flip(e)) === tail(e)
        - !is_primary(e) || (!isnothing(tail(e)) && !isnothing(head(e)) && (!isnothing(left(e)) || !isnothing(right(e))))
        - !is_dual(e) || (!isnothing(left(e)) && !isnothing(right(e)) && (!isnothing(head(e)) || !isnothing(tail(e))))
        """
        function _check_consistency(mesh::PolyhedralMesh)
            for (i, e) in enumerate(vcat(collect(values(edges(mesh))), collect(values(dual_edges(mesh)))))
                if !(tail(next(e)) === tail(e))
                    @warn "e: $e"
                    @warn "tail(next(e)) != tail(e): $(tail(next(e))) !== $(tail(e))"
                    return false
                elseif !(left(e) === right(next(e)))
                    @warn "e: $e"
                    @warn "left(e) != right(next(e)): $(left(e)) !== $(right(next(e)))"
                    return false
                elseif !(tail(flip(e)) === tail(e))
                    @warn "e: $e"
                    @warn "tail(flip(e)) != tail(e): $(tail(flip(e))) !== $(tail(e))"
                    return false
                end
            end
            return true
        end

        tetra = tetrahedron()
        @test _check_consistency(tetra)
        tetra_with_hole = PolyhedralMesh([[1, 2, 3], [1, 2, 4], [1, 3, 4]])
        @test _check_consistency(tetra_with_hole)
        @test_throws ArgumentError PolyhedralMesh([[1, 2, 3], [1, 2, 4], [1, 2, 5]])
        octa = octahedron()
        @test _check_consistency(octa)

        for f in faces(tetra)
            @test length(vertices(f)) == 3
        end

        for f in faces(octa)
            @test length(vertices(f)) == 3
        end
    end
end

@testset "Rigidity Tests" begin
    @testset "Kabsch Algorithm: 100 Points in R³ (Perfect Fit)" begin
        ## 1. SETUP: Define problem size and parameters
        num_points = 100
        dimension = 3

        ## 2. GENERATE a known "true" transformation
        # Create a random, proper rotation matrix
        R_true, _ = qr(randn(dimension, dimension))
        if det(R_true) < 0
            # Ensure it's a proper rotation (determinant of +1)
            R_true = R_true * Diagonal([1, 1, -1])
        end
        # Create a random translation vector
        t_true = 100 .* (rand(dimension) .- 0.5)

        ## 3. CREATE the initial set of "moving" points (P)
        P = 20 .* randn(dimension, num_points)

        ## 4. CREATE the "reference" point set (Q)
        # Apply the exact transformation with no noise.
        Q = R_true * P .+ t_true

        ## 5. RUN the Kabsch algorithm
        (R_calc, t_calc, rmsd_val) = kabsch(P, Q)

        ## 6. VERIFY the results
        # The calculated transformation should be nearly identical to the true one.
        @test R_calc ≈ R_true
        @test t_calc ≈ t_true

        # The RMSD should be effectively zero since there is no noise.
        @test rmsd_val ≈ 0.0 atol = 1e-12
    end
end